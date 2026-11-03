"""
参数扫描脚本：对 比例([B]/max_B : TDF)、max_B、C 链 Kd_C 三维批量扫描，汇总最终产物分布。

扫描维度：
  1. ratio（每点 B 与 TDF 的比 = [B]/max_B : TDF），con_B_per_site = ratio * con_TDF
  2. max_B（主链最大 B 数），con_B = max_B * con_B_per_site
  3. Kd_C（C 位点解离常数），从 20 uM 到 50 nM 对数取点

输出（全部写入 source_data/ 子目录；文件名带比例标签与 max_B 条件）：
  <ratio>_scan_summary.csv       —— 每个 (比例, max_B, Kd_C) 一行，关键汇总指标
  <ratio>_scan_final_states.csv  —— 长格式，每个组合的最终物种分布
  <ratio>_maxB9_n9_TDF_yield.csv —— max_B=9 时 TDF 各状态占比 vs Kd_C
  <ratio>_maxB9_n9_TDF_yield.png —— 对应 stacked 面积图
  kinetics_<ratio>_maxB{n}_KdC{m}nM.csv —— 每个组合的完整时间序列（每个时刻的状态）

用法：
  python scan.py           # 默认 t_max=10000 s
  python scan.py 1000      # 覆盖 t_max
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tDFTM_network_SSA as m

# ==================== 扫描参数 ====================
con_TDF = m.con_TDF          # TDF 初始浓度 (M)，默认 100 nM
# 比例：[B]/max_B : TDF（每点 B 与 TDF 的比）。con_B_per_site = ratio * con_TDF
ratio_list = [0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
              2.0, 2.25, 2.5, 2.75, 3.0, 5.0, 7.5, 10.0]
max_B_list = [3]             # 扫描的主链最大 B 数
Kd_C_list = [5e-5, 2e-5, 1e-5, 5e-6, 2e-6, 1e-6, 5e-7, 2e-7, 1e-7, 5e-8]  # 50 uM .. 50 nM
t_max = 10000.0              # 模拟时长（足够 C 链接近平衡）
SEED = None                  # None = 每次随机；设整数可复现
OUT_DIR = 'source_data_1030' # 输出子目录
# ==================================================

N_A = m.N_A
V = 1.0e-11                  # 模拟盒子体积 (L)；只影响噪声尺度与运行耗时，不影响平均值
Kon_B = m.Kon_B
Koff_B = m.Koff_B
Kon_C = m.Kon_C

k_fB = Kon_B / (N_A * V)
k_rB = Koff_B
k_fC = Kon_C / (N_A * V)


def ratio_label(r):
    """把比例 r = ([B]/max_B)/con_TDF 转成文件名标签（TDF : 每点B = 1 : r）"""
    return f'TDF1_B{r:g}'


def make_n9_yield(final_rows, out_dir=OUT_DIR, rl=''):
    """从长格式最终分布提取 max_B=9 的 TDF yield，输出 CSV 与 stacked 图"""
    os.makedirs(out_dir, exist_ok=True)
    prefix = f'{rl}_' if rl else ''
    df = pd.DataFrame(final_rows)
    n9 = df[df['max_B'] == 9]
    if n9.empty:
        print('无 max_B=9 数据，跳过 yield 汇总')
        return

    piv = n9.pivot_table(index='Kd_C_nM', columns='species',
                         values='count', aggfunc='sum', fill_value=0)
    non_B = [c for c in piv.columns if c != 'B']
    total_TDF = piv[non_B].sum(axis=1)
    frac = piv[non_B].div(total_TDF, axis=0)

    keep = ['TDFB9'] + [f'TDFB9C{j}' for j in range(1, m.j_max(9) + 1)]
    keep = [c for c in keep if c in frac.columns]
    y = frac[keep].sort_index(ascending=False)

    y.to_csv(os.path.join(out_dir, f'{prefix}maxB9_n9_TDF_yield.csv'))

    fig, ax = plt.subplots(figsize=(10, 6))
    x = y.index.values.astype(float)
    ax.stackplot(x, y.values.T, labels=keep,
                 colors=plt.cm.viridis(np.linspace(0, 1, len(keep))), alpha=0.85)
    ax.set_xscale('log')
    ax.invert_xaxis()
    ax.set_xlabel('Kd_C (nM)')
    ax.set_ylabel('Fraction of TDF')
    ax.set_title(f'TDF-B9 yield distribution vs Kd_C ({rl})')
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    png_path = os.path.join(out_dir, f'{prefix}maxB9_n9_TDF_yield.png')
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    print(f'n9_TDF_yield 已保存: {os.path.join(out_dir, f"{prefix}maxB9_n9_TDF_yield.csv")} 和 {png_path}')


def main(t_max=t_max):
    os.makedirs(OUT_DIR, exist_ok=True)

    for ratio in ratio_list:
        con_B_per_site = ratio * con_TDF
        rl = ratio_label(ratio)
        summary_rows = []
        final_rows = []

        for max_B in max_B_list:
            con_B = max_B * con_B_per_site
            names, index = m.build_species(max_B)

            n_TDF = int(round(con_TDF * N_A * V))
            n_B = int(round(con_B * N_A * V))
            init = np.zeros(len(names), dtype=np.int64)
            init[index['TDF']] = n_TDF
            init[index['B']] = n_B

            for Kd_C in Kd_C_list:
                Koff_C = Kd_C * Kon_C
                reactions = m.build_reactions(index, max_B, k_fB, k_rB, k_fC, Koff_C)

                print(f'\n===== 条件: ratio={rl}, max_B={max_B}, Kd_C={Kd_C*1e6:.2f} uM, con_B={con_B*1e9:.0f} nM =====', flush=True)

                result, n_events, t_end = m.gillespie(reactions, init, t_max, seed=SEED, progress_dt=100)
                final = result[-1, 1:].astype(np.int64)

                # 保存完整时间序列（每个时刻的状态），与原始程序 KMC_TDF.py 输出格式一致：
                # 每行 = 一个时刻 [time, TDF, B, TDFB1, ... , TDFB9C7]
                Kd_nM = int(round(Kd_C * 1e9))
                ts_path = os.path.join(OUT_DIR, f'kinetics_{rl}_maxB{max_B}_KdC{Kd_nM}nM.csv')
                pd.DataFrame(result, columns=['time'] + names).to_csv(ts_path, index=False)

                primary_B = 0
                C_bound = 0
                with_C = 0
                for i in range(max_B + 1):
                    for j in range(m.j_max(i) + 1):
                        cnt = int(final[index[m.species_name(i, j)]])
                        primary_B += i * cnt
                        C_bound += j * cnt
                        if j >= 1:
                            with_C += cnt

                avg_primary_B = primary_B / n_TDF
                avg_C = C_bound / n_TDF
                frac_with_C = with_C / n_TDF
                bare_full = int(final[index[m.species_name(max_B, 0)]]) / n_TDF
                free_B = int(final[index['B']])

                summary_rows.append({
                    'ratio': ratio,
                    'B_per_site_nM': con_B_per_site * 1e9,
                    'max_B': max_B,
                    'con_B_nM': con_B * 1e9,
                    'Kd_C_nM': Kd_C * 1e9,
                    'avg_primary_B': avg_primary_B,
                    'avg_C_per_TDF': avg_C,
                    'frac_with_C': frac_with_C,
                    'frac_bare_TDFBmax': bare_full,
                    'free_B': free_B,
                    'n_events': n_events,
                    't_end': t_end,
                })

                for name in names:
                    cnt = int(final[index[name]])
                    if cnt > 0:
                        final_rows.append({
                            'ratio': ratio,
                            'max_B': max_B,
                            'con_B_nM': con_B * 1e9,
                            'Kd_C_nM': Kd_C * 1e9,
                            'species': name,
                            'count': cnt,
                        })

                print(f'ratio={rl}, max_B={max_B}, Kd_C={Kd_C*1e6:6.2f} uM -> '
                      f'avg_C={avg_C:.3f}, frac_with_C={frac_with_C:.3f}')

        summary_name = f'{rl}_scan_summary.csv'
        final_name = f'{rl}_scan_final_states.csv'
        pd.DataFrame(summary_rows).to_csv(os.path.join(OUT_DIR, summary_name), index=False)
        pd.DataFrame(final_rows).to_csv(os.path.join(OUT_DIR, final_name), index=False)
        print(f'\n汇总已保存: {os.path.join(OUT_DIR, summary_name)}')
        print(f'分布已保存: {os.path.join(OUT_DIR, final_name)}')

        make_n9_yield(final_rows, OUT_DIR, rl)


if __name__ == '__main__':
    t = float(sys.argv[1]) if len(sys.argv) > 1 else t_max
    main(t)