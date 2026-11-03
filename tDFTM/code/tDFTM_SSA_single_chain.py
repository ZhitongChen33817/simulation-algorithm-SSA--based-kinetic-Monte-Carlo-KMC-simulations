import math
import random
import numpy as np
from typing import List
import pandas as pd

def read_order(path,volume_so):
    df2 = pd.read_csv(path, sep=',', header=None, usecols=[0], names=[0])
    Read_forward = df2[0].values
    print(f"component: {len(Read_forward)}, Volume: {volume_so / 1e-6:.1e} uL,"
          f"K_f {Read_forward[0]:.1e}, K_r: {Read_forward[1]:.1e}, "
          f"K_d: {Read_forward[1]/(Read_forward[0]*volume_so):.1e}")
    return  Read_forward


def generate_reactions(seq,cnt):
    """生成所有9个反应方向（18对可逆反应）"""
    species = {
        'A': 0, 'B': 1,'AB': 2,  'AB2': 3,  'AB3': 4,  'AB4': 5,   'AB5': 6,
        'AB6': 7,  'AB7': 8,  'AB8': 9,   'AB9': 10   }

    reactions = [
        # 1. A + B ⇌ AB
        ("A+B→AB", seq[0] * cnt[0] * cnt[1], [-1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
        ("AB→A+B", seq[1] * cnt[2], [1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0]),
        # 2. AB + B ⇌ AB2
        ("AB+B→AB2", seq[2] * cnt[2] * cnt[1], [0, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0]),
        ("AB2→AB+B", seq[3] * cnt[3], [0, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0]),   
        # 3. AB2 + B ⇌ AB3
        ("AB2+B→AB3", seq[4] * cnt[3] * cnt[1], [0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0]),
        ("AB3→AB2+B", seq[5] * cnt[4], [0, 1, 0, 1, -1, 0, 0, 0, 0, 0, 0]),   
        # 4. AB3 + B ⇌ AB4
        ("AB3+B→AB4", seq[6] * cnt[4] * cnt[1], [0, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0]),
        ("AB4→AB3+B", seq[7] * cnt[5], [0, 1, 0, 0, 1, -1, 0, 0, 0, 0, 0]),   
        # 5. AB4 + B ⇌ AB5
        ("AB4+B→AB5", seq[8] * cnt[5] * cnt[1], [0, -1, 0, 0, 0, -1, 1, 0, 0, 0, 0]),
        ("AB5→AB4+B", seq[9] * cnt[6], [0, 1, 0, 0, 0, 1, -1, 0, 0, 0, 0]),   
        # 6. AB5 + B ⇌ AB6
        ("AB5+B→AB6", seq[10] * cnt[6] * cnt[1], [0, -1, 0, 0, 0, 0, -1, 1, 0, 0, 0]),
        ("AB6→AB5+B", seq[11] * cnt[7], [0, 1, 0, 0, 0, 0, 1, -1, 0, 0, 0]),   
        # 7. AB6 + B ⇌ AB7
        ("AB6+B→AB7", seq[12] * cnt[7] * cnt[1], [0, -1, 0, 0, 0, 0, 0, -1, 1, 0, 0]),
        ("AB7→AB6+B", seq[13] * cnt[8], [0, 1, 0, 0, 0, 0, 0, 1, -1, 0, 0]),   
        # 8. AB7 + B ⇌ AB8
        ("AB7+B→AB8", seq[14] * cnt[8] * cnt[1], [0, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0]),
        ("AB8→AB7+B", seq[15] * cnt[9], [0, 1, 0, 0, 0, 0, 0, 0, 1, -1, 0]),
        # 9. AB8 + B ⇌ AB9
        ("AB8+B→AB9", seq[16] * cnt[9] * cnt[1], [0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 1]),
        ("AB9→AB8+B", seq[17] * cnt[10], [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, -1])
    ]
    return reactions


def gillespie_simulation(seq_k2: List[float], initial_counts: List[int],t_max: float):
    """
    完整版Gillespie模拟
    参数:
        c: 132个反应参数 [c1_forward, c1_backward, ..., c66_forward, c66_backward]
        initial_counts: 初始分子数 (长度27的列表)
        t_max: 最大模拟时间
    """

    result_output = []
    results_list = []
    counts = initial_counts.copy()
    t = 0.0
    seq_reaction = 0
    next_save_time = 0.025
    reactions = generate_reactions(seq_k2, initial_count)

    # 首次输出
    results_list.append(np.hstack([t, counts]))
    print(f"start: Time: {t:.1f}, A: {(counts[0] / 6.02e23 / 1e-11):.1e}, B: {(counts[1] / 6.02e23 / 1e-11):.1e}")

    while t < t_max:
        # 计算所有反应倾向
        reaction_1 = generate_reactions(seq_k2, counts)
        sum_column2 = sum(row[1] for row in reaction_1)

        if sum_column2 <= 1e-12:  # 防止浮点误差
            print("所有反应倾向为零，模拟终止")
            print('seq:', seq_reaction / 1e3, 't:', t)
            results_list.append(np.hstack([t, counts]))
            result_output = np.array(results_list)
            break

        # 生成下一事件时间
        tau = -math.log(random.random()) / sum_column2
        t += tau
        seq_reaction += 1
        # 处理输出

        if t >= next_save_time:
            print('seq:', seq_reaction / 5e4, 't:', t)
            results_list.append(np.hstack([t, counts]))
            result_output = np.array(results_list)

            # 根据当前时间动态调整下次记录间隔
            if t < 30:
                next_save_time += 0.025  # 30秒内每0.1秒记录
            elif t < 200:
                next_save_time += 1  # 30-200秒每1秒记录
            elif t < 2000:
                next_save_time += 10  # 200-2000秒每1秒记录

            else:
                next_save_time += 100  # 2000秒以上每10秒记录

        # 选择反应
        r = random.random() * sum_column2
        cum_sum = 0.0
        for i in range(len(reactions)):
            cum_sum += reaction_1[i][1]
            if cum_sum >= r:
                # 更新分子数
                delta = reactions[i][2]
                # 先计算新值但不立即更新
                new_counts = [x + y for x, y in zip(counts, delta)]
                if all(c >= 0 for c in new_counts):  # 注意这里用all，因为要求所有都非负
                    counts = new_counts
                break
    return result_output


#Kon= 1.66e6;  Kd = 10 nM(Kd= Koff/Kon); Koff= 1.66e-2
#一级反应速率的hc_f= Koff/N_A(2.74e-26), 二级反应速率的hc_r = Kon/(N_A*V) (2.76e-18 / 1e-11 = 2.76e-7),
N_A = 6.022e23
n = 8
con_A = 100e-9 *N_A  #100 nM
con_B = n * 150e-9 *N_A  #100 nM
Volume = 1e-11#1uL 体系
initial_count = np.array([con_A*Volume, con_B*Volume] + [0] * 9) #A与B分子的初始分子数； 100 nM *1 uL*NA = 6e11

order_path = 'input/seq_tDF_9unit.csv'
seq_k =read_order(order_path,Volume)
print("seq:", len(seq_k))


result_TDF = gillespie_simulation(seq_k, initial_count,1000)
filename = 'onehalf_KMC_TDF_100nM_data_8unit.csv'  # 定义文件名变量
np.savetxt(filename, result_TDF, delimiter=',')
print(f"数据已保存至文件: {filename}")  # 打印文件名
