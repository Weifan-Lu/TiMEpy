import numpy as np
import matplotlib.pyplot as plt


def acoth(x):
    """
    计算反双曲余切：acoth(x) = 0.5 * log((x+1)/(x-1))
    注意：仅对 |x| > 1 定义。
    """
    return 0.5 * np.log((x + 1) / (x - 1))


def b_positive_discrete(magnitudes, delta=0.05):
    """
    利用 b-positive 方法估计 b 值（针对离散化震级数据）。

    参数:
        magnitudes (array-like): 按时间顺序排列的地震震级序列。
        delta (float): 离散化参数（震级以 2*delta 为间隔离散，默认 delta=0.05 对应步长 0.1）。

    返回:
        float: 估计的 b 值，如果计算失败则返回 None。

    逻辑说明:
        1. 首先将震级数据离散化：按照 2*delta 的步长进行四舍五入。
           例如，delta=0.05 时，2*delta=0.1，即保留 1 位小数。
        2. 计算相邻震级之间的差值。
        3. 设 Mc' = 2*delta，筛选出差值大于或等于 Mc'（考虑数值容差）的数据。
        4. 计算筛选后差值的均值，利用公式：

            b = [1/(delta * ln(10))] * acoth( (mean(diff) - Mc' + delta)/delta )

           其中，acoth 为反双曲余切函数。
    """
    # 根据 2*delta 的数量级确定需要保留的小数位数
    # 例如：delta=0.05 -> 2*delta=0.1 -> -log10(0.1)=1
    decimals = int(round(-np.log10(2 * delta)))
    discretized = np.round(magnitudes, decimals=decimals)

    # 计算相邻震级差值
    diffs = np.diff(discretized)

    # 设置 Mc' = 2*delta，并使用一个小的容差
    Mc_prime = 2 * delta
    tol = 1e-6
    # 筛选出大于或等于 Mc' 的差值
    pos_diffs = diffs[diffs >= Mc_prime - tol]

    if len(pos_diffs) == 0:
        print("没有找到满足条件的正震级差值 (>= Mc').")
        b_value = 1000
        return None

    mean_diff = np.mean(pos_diffs)
    # 计算 acoth 的输入参数：(mean_diff - Mc' + delta)/delta
    arg = (mean_diff - Mc_prime + delta) / delta
    if arg <= 1:
        print("acoth 的参数小于等于 1，无法计算 acoth。")
        b_value = 1000
        return None

    # 根据 MATLAB 公式计算 b 值
    b_value = 1.0 / (delta * np.log(10)) * acoth(arg)
    return b_value


def b_positive_discrete_bootstrap(magnitudes, delta=0.05, n_boot=1000):
    """
    利用 bootstrap 方法对 b_positive_discrete 的估计误差进行估计。

    参数:
        magnitudes (array-like): 原始震级序列（按时间顺序排列）。
        delta (float): 离散化参数。
        n_boot (int): bootstrap 重采样次数。

    返回:
        tuple: (b_mean, b_std)
            b_mean - bootstrap 平均估计的 b 值
            b_std  - b 值的标准差，作为误差估计
    """
    boot_b_values = []
    n = len(magnitudes)

    for i in range(n_boot):
        # 对原始目录进行有放回重采样，并保持时间顺序
        indices = np.random.choice(np.arange(n), size=n, replace=True)
        indices.sort()  # 保证重采样后数据依然按时间排序
        sample = magnitudes[indices]
        b_val_sample = b_positive_discrete(sample, delta=delta)
        if b_val_sample is not None:
            boot_b_values.append(b_val_sample)

    boot_b_values = np.array(boot_b_values)
    if len(boot_b_values) == 0:
        print("Bootstrap 未能获得有效的 b 值估计。")
        return None,None

    return np.mean(boot_b_values), np.std(boot_b_values)

#
# # Data processing
# filename = '/Users/luwf/PycharmProjects/UTokyoEPS/work_file/Tidal_Triggering/TME_package/ex_noto/output/catalog/New_Noto_select_catalog.txt'
# data = np.loadtxt(filename)
#
# # 按时间排序（假设数据的前6列依次为 年, 月, 日, 时, 分, 秒）
# sort_indices = np.lexsort((data[:, 5], data[:, 4], data[:, 3],
#                            data[:, 2], data[:, 1], data[:, 0]))
# sorted_data = data[sort_indices]
#
# # 提取震级数据（假定震级位于第10列，即索引9）
# magnitudes = sorted_data[:, 9]
#
# # 计算 b-positive 值（离散化方法）
# b_val = b_positive_discrete(magnitudes, delta=0.05)
#
# if b_val is not None:
#     print(f"Estimated b-value: {b_val:.2f}")
# else:
#     print("Failed to calculate b-value.")
#
# # 利用 bootstrap 方法估计 b 值的误差
# b_mean, b_std = b_positive_discrete_bootstrap(magnitudes, delta=0.05, n_boot=1000)
# if b_mean is not None:
#     print(f"Bootstrap estimated b-value: {b_mean:.2f} ± {b_std:.2f}")
# else:
#     print("Bootstrap estimation failed.")


# 可选：绘制 b 值分布直方图（基于 bootstrap 样本）
# 若希望可视化 bootstrap 得到的 b 值分布，可将 bootstrap 的 b 值数组保留并绘图
# 下面给出一个简单示例：
def get_bootstrap_b_values(magnitudes, delta=0.05, n_boot=1000):
    boot_b_values = []
    n = len(magnitudes)
    for i in range(n_boot):
        indices = np.random.choice(np.arange(n), size=n, replace=True)
        indices.sort()
        sample = magnitudes[indices]
        b_val_sample = b_positive_discrete(sample, delta=delta)
        if b_val_sample is not None:
            boot_b_values.append(b_val_sample)
    return np.array(boot_b_values)
