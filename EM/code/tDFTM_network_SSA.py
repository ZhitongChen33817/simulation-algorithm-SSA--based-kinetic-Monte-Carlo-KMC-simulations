
import math
import sys
import numpy as np
import pandas as pd

# ==================== 可配置参数 ====================
N_A = 6.022e23          # 阿伏伽德罗常数 (1/mol)
V = 1.0e-11            # 体积 (L)；生产运行设置，只影响噪声尺度与速度，不影响平均结果
max_B = 9               # 主链最大 B 数（可截断到任意 <=8）

# B 位点（主链）参数（Kd = 20 pM）
Kon_B = 1.66e6          # M^-1 s^-1
Koff_B = 3.32e-5        # s^-1

# C 位点（次链）参数（Kd = 0.1 uM）
Kon_C = 5.0e2           # M^-1 s^-1
Koff_C = 5.0e-5         # s^-1

# 归一化（标准约定）
k_fB = Kon_B / (N_A * V)   # 主链正反应（二级）
k_rB = Koff_B              # 主链逆反应（一级）
k_fC = Kon_C / (N_A * V)   # 次链正反应（二级）
k_rC = Koff_C              # 次链逆反应（一级）

# 初始浓度 (M)
con_TDF = 100e-9          # 100 nM
con_B = 9 * 300e-9        # 1200 nM（12 倍过量；每个 TDF 最多结合 14 个 B）

# 模拟时长 (s)，可用命令行参数覆盖：python tDFTM_network_SSA.py 10
t_max = 1000.0

# 随机种子（None = 每次不同；设整数可复现）
SEED = None
# ==================================================


def j_max(i):
    """主链 B 数为 i 时，次链 C 的最大数（C 位点逻辑：j <= i-2）"""
    return max(0, i - 2)


def species_name(i, j):
    """把 (i, j) 映射为物种名"""
    if i == 0 and j == 0:
        return 'TDF'
    if j == 0:
        return f'TDFB{i}'
    return f'TDFB{i}C{j}'


def build_species(max_B):
    """生成物种名列表与索引字典"""
    names = ['TDF', 'B']
    index = {'TDF': 0, 'B': 1}
    # 主链（j = 0）
    for i in range(1, max_B + 1):
        name = f'TDFB{i}'
        index[name] = len(names)
        names.append(name)
    # 次链（j >= 1）
    for i in range(3, max_B + 1):
        for j in range(1, i - 1):   # j = 1 .. i-2
            name = f'TDFB{i}C{j}'
            index[name] = len(names)
            names.append(name)
    return names, index


def build_reactions(index, max_B, k_fB, k_rB, k_fC, k_rC):
    """生成所有反应（向量化表示：k, r1, r2, delta）"""
    n = len(index)
    rxn_names = []
    k_list = []
    r1_list = []
    r2_list = []
    delta_list = []
    bi = index['B']

    def add(name, k, order, reactants, changes):
        rxn_names.append(name)
        k_list.append(k)
        r1_list.append(reactants[0])
        r2_list.append(reactants[1] if order == 2 else -1)
        d = np.zeros(n, dtype=np.int64)
        for idx, ch in changes:
            d[idx] += ch
        delta_list.append(d)

    # 主链：TDFB{i-1}C{j} + B <-> TDFB{i}C{j}
    for i in range(1, max_B + 1):
        for j in range(j_max(i - 1) + 1):
            src = species_name(i - 1, j)
            dst = species_name(i, j)
            si = index[src]
            di = index[dst]
            add(f'{src} + B -> {dst}', k_fB, 2, (si, bi), [(si, -1), (bi, -1), (di, +1)])
            add(f'{dst} -> {src} + B', k_rB, 1, (di,), [(di, -1), (si, +1), (bi, +1)])

    # 次链：TDFB{i}C{j-1} + B <-> TDFB{i}C{j}
    for i in range(3, max_B + 1):
        for j in range(1, i - 1):   # j = 1 .. i-2
            src = species_name(i, j - 1)
            dst = species_name(i, j)
            si = index[src]
            di = index[dst]
            add(f'{src} + B -> {dst}', k_fC, 2, (si, bi), [(si, -1), (bi, -1), (di, +1)])
            add(f'{dst} -> {src} + B', k_rC, 1, (di,), [(di, -1), (si, +1), (bi, +1)])

    k = np.array(k_list)
    r1 = np.array(r1_list)
    r2 = np.array(r2_list)
    delta = np.array(delta_list)
    return rxn_names, k, r1, r2, delta


def gillespie(reactions, initial_counts, t_max, seed=None, progress_dt=None):
    """通用 Gillespie SSA（直接法）。progress_dt 给定则每 progress_dt 秒打印一次进度。"""
    rxn_names, k, r1, r2, delta = reactions
    if seed is not None:
        np.random.seed(seed)

    counts = initial_counts.astype(np.int64).copy()
    t = 0.0
    n_events = 0

    results = [np.concatenate([[t], counts])]

    def save_interval(t):
        if t < 30:
            return 0.025
        elif t < 200:
            return 1.0
        elif t < 2000:
            return 10.0
        else:
            return 100.0

    next_save = 0.025
    next_progress = progress_dt if progress_dt is not None else None

    while t < t_max:
        # 向量化倾向性
        a = k * counts[r1]
        second = r2 >= 0
        a[second] *= counts[r2[second]]
        a0 = a.sum()

        if a0 <= 0.0:
            break

        tau = -math.log(np.random.random()) / a0
        t += tau

        r = np.random.random() * a0
        cum = np.cumsum(a)
        idx = min(int(np.searchsorted(cum, r)), len(a) - 1)
        counts += delta[idx]
        n_events += 1

        if t >= next_save:
            results.append(np.concatenate([[t], counts]))
            next_save += save_interval(next_save)

        if next_progress is not None and t >= next_progress:
            print(f'  进度: t={t:.0f} s, 事件数={n_events}', flush=True)
            next_progress += progress_dt

    results.append(np.concatenate([[t], counts]))
    return np.array(results), n_events, t


def main(t_max=t_max):
    print('=== TDF 双位点顺序结合模型 ===')
    names, index = build_species(max_B)
    reactions = build_reactions(index, max_B, k_fB, k_rB, k_fC, k_rC)
    rxn_names, k, r1, r2, delta = reactions

    print(f'物种数: {len(names)}，单向反应数: {len(rxn_names)}')
    print(f'归一化速率常数: k_fB={k_fB:.3e}, k_rB={k_rB:.3e}, k_fC={k_fC:.3e}, k_rC={k_rC:.3e}')

    # 初始计数
    n_TDF = int(round(con_TDF * N_A * V))
    n_B = int(round(con_B * N_A * V))
    initial_counts = np.zeros(len(names), dtype=np.int64)
    initial_counts[index['TDF']] = n_TDF
    initial_counts[index['B']] = n_B
    print(f'初始: TDF={n_TDF}, B={n_B}')

    # 跑模拟
    result, n_events, t_end = gillespie(reactions, initial_counts, t_max, seed=SEED)
    print(f'模拟结束: t={t_end:.2f} s, 事件数={n_events}')

    # 守恒检查
    final_state = result[-1, 1:].astype(np.int64)
    tdf_total = int(final_state.sum() - final_state[index['B']])
    print(f'TDF 守恒: 初始 {n_TDF} -> 最终 {tdf_total}  ({ "一致" if tdf_total == n_TDF else "不一致!" })')

    # 保存 CSV
    header = ['time'] + names
    out = 'tDFTM_network_result.csv'
    pd.DataFrame(result, columns=header).to_csv(out, index=False)
    print(f'结果已保存: {out}')

    # 最终分布（按数量降序，前 12）
    print('\n最终分布（按数量降序，前 12）:')
    order = np.argsort(-final_state)
    for idx in order[:12]:
        print(f'  {names[idx]:<10} {int(final_state[idx])}')


if __name__ == '__main__':
    t = float(sys.argv[1]) if len(sys.argv) > 1 else t_max
    main(t)
