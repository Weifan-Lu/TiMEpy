import numpy as np
from scipy.optimize import minimize
import visualization as vis
import matplotlib.pyplot as plt


def load_txt_data(filename):
    """
    读取 txt 文件，跳过第一行 header，
    文件采用制表符分隔，返回二维 NumPy 数组
    """
    return np.loadtxt(filename, delimiter='\t', skiprows=1)

#
#
# def calculate_tidal_phase(ttide_cut, Stress_cut, t_ref):
#     """
#     根据给定的时间窗口 ttide_cut 和应力序列 StressS_cut，
#     以及地震事件时间 t_ref（均为 datenum 格式），计算潮汐相位等信息。
#
#     输入:
#       ttide_cut   : 1D numpy 数组，时间序列（单位：天）
#       StressS_cut : 1D numpy 数组，对应的应力序列
#       t_ref       : 地震事件时间（天，datenum 格式）
#     输出:
#       ephase           : 潮汐相位（单位：弧度）
#       maxtime          : 窗口内局部极大值时刻
#       mintime1         : 左侧局部极小时刻
#       mintime2         : 右侧局部极小时刻
#       stress_at_t      : 地震时刻对应的应力值
#       stress_rate_at_t : 地震时刻对应的应力变化率
#     """
#     # 计算应力变化率（差分除以时间间隔）
#     stress_rate = np.diff(Stress_cut) / np.diff(ttide_cut)
#
#     # 找到离 t_ref 最近的索引
#     i = np.argmin(np.abs(ttide_cut - t_ref))
#     stress_at_t = Stress_cut[i]
#     if i < len(stress_rate):
#         stress_rate_at_t = stress_rate[i]
#     else:
#         stress_rate_at_t = np.nan
#
#     ephase = None
#     maxtime = None
#     mintime1 = None
#     mintime2 = None
#
#     if stress_rate[i] >= 0:
#         # 正应力变化率情况
#         j = i
#         while j < len(stress_rate) and stress_rate[j] >= 0:
#             j += 1
#         maxtime = ttide_cut[j - 1] if j - 1 >= 0 else ttide_cut[i]
#
#         k = j
#         while k < len(stress_rate) and stress_rate[k] < 0:
#             k += 1
#         mintime2 = ttide_cut[k - 1] if k - 1 >= 0 else ttide_cut[i]
#
#         l = i
#         while l > 0 and stress_rate[l] >= 0:
#             l -= 1
#         mintime1 = ttide_cut[l + 1] if l + 1 < len(ttide_cut) else ttide_cut[i]
#
#         # 计算潮汐相位：遍历 m 从 0 到 180（步长 0.5°）
#         for m in np.arange(0, 180.5, 0.5):
#             ephasetime = (maxtime - mintime1) / 180.0 * m + mintime1
#             if ttide_cut[i] <= ephasetime:
#                 ephase = -180 + m
#                 break
#     else:
#         # 负应力变化率情况
#         j = i
#         while j < len(stress_rate) and stress_rate[j] < 0:
#             j += 1
#         mintime2 = ttide_cut[j - 1] if j - 1 >= 0 else ttide_cut[i]
#
#         k = i
#         while k > 0 and stress_rate[k] < 0:
#             k -= 1
#         maxtime = ttide_cut[k + 1] if k + 1 < len(ttide_cut) else ttide_cut[i]
#
#         l = k
#         while l > 0 and stress_rate[l] >= 0:
#             l -= 1
#         mintime1 = ttide_cut[l + 1] if l + 1 < len(ttide_cut) else ttide_cut[i]
#
#         for m in np.arange(0, 180.5, 0.5):
#             ephasetime = (mintime2 - maxtime) / 180.0 * m + maxtime
#             if ttide_cut[i] <= ephasetime:
#                 ephase = m
#                 break
#
#     # 转换为弧度
#     if ephase is not None:
#         ephase = np.deg2rad(ephase)
#     else:
#         ephase = np.nan
#
#     return ephase, maxtime, mintime1, mintime2, stress_at_t, stress_rate_at_t


def calculate_tidal_phase(ttide_cut, StressS_cut, t_ref):
    """
    Calculate tidal phase and related information.

    Inputs:
      ttide_cut   : 1D numpy array, time series
      StressS_cut : 1D numpy array, tidal stress sequence
      t_ref       : Earthquake event occurrence time
    Outputs:
      ephase           : Tidal phase (in radians)
      maxtime          : Time corresponding to the local maximum
      mintime1         : Time corresponding to the local minimum on the left
      mintime2         : Time corresponding to the local minimum on the right
      stress_at_t      : Tidal stress value at the earthquake time
      stress_rate_at_t : Tidal stress rate at the earthquake time
    """
    # Step 1: Compute the stress rate
    stress_rate = np.diff(StressS_cut) / np.diff(ttide_cut)

    # Step 2: Find the index of the earthquake event time (closest index)
    i = np.argmin(np.abs(ttide_cut - t_ref))

    # Step 3: Get the tidal stress value and stress rate at the earthquake time
    stress_at_t = StressS_cut[i]
    if i < len(stress_rate):
        stress_rate_at_t = stress_rate[i]
    else:
        stress_rate_at_t = np.nan

    # Initialize output variables
    ephase = None
    maxtime = None
    mintime1 = None
    mintime2 = None

    # Step 4: Determine the trend of the stress rate
    if stress_rate[i] >= 0:
        # Positive stress rate case
        # Find the local maximum
        j = i
        while j < len(stress_rate) and stress_rate[j] >= 0:
            j += 1
        maxtime = ttide_cut[j - 1]  # Time corresponding to the local maximum

        # Find the local minimum on the right
        k = j
        while k < len(stress_rate) and stress_rate[k] < 0:
            k += 1
        mintime2 = ttide_cut[k - 1]  # Time corresponding to the local minimum on the right

        # Find the local minimum on the left
        l = i
        while l > 0 and stress_rate[l] >= 0:
            l -= 1
        mintime1 = ttide_cut[l + 1]  # Time corresponding to the local minimum on the left

        # Compute the tidal phase
        for m in np.arange(0, 180.5, 0.5):
            ephasetime = (maxtime - mintime1) / 180.0 * m + mintime1
            if ttide_cut[i] <= ephasetime:
                ephase = -180 + m
                break
    else:
        # Negative stress rate case
        # Find the local minimum on the right
        j = i
        while j < len(stress_rate) and stress_rate[j] < 0:
            j += 1
        mintime2 = ttide_cut[j - 1]  # Time corresponding to the local minimum on the right

        # Find the local maximum on the left
        k = i
        while k > 0 and stress_rate[k] < 0:
            k -= 1
        maxtime = ttide_cut[k + 1]  # Time corresponding to the local maximum

        # Find the local minimum on the left
        l = k
        while l > 0 and stress_rate[l] >= 0:
            l -= 1
        mintime1 = ttide_cut[l + 1]  # Time corresponding to the local minimum on the left

        # Compute the tidal phase
        for m in np.arange(0, 180.5, 0.5):
            ephasetime = (mintime2 - maxtime) / 180.0 * m + maxtime
            if ttide_cut[i] <= ephasetime:
                ephase = m
                break

    # Convert tidal phase from degrees to radians
    if ephase is not None:
        ephase = np.deg2rad(ephase)
    else:
        ephase = np.nan

    return ephase, maxtime, mintime1, mintime2, stress_at_t, stress_rate_at_t


def calculate_fast_tidal_phase(ttide_cut, Stress_cut, t_ref):
    """
    根据给定的时间窗口 ttide_cut 和应力序列 Stress_cut，
    以及地震事件时间 t_ref（均为 datenum 格式），计算潮汐相位等信息。

    输入:
      ttide_cut   : 1D numpy 数组，时间序列（单位：天）
      Stress_cut  : 1D numpy 数组，对应的应力序列
      t_ref       : 地震事件时间（天，datenum 格式）
    输出:
      ephase           : 潮汐相位（单位：弧度）
      maxtime          : 窗口内局部极大值时刻
      mintime1         : 左侧局部极小时刻
      mintime2         : 右侧局部极小时刻
      stress_at_t      : 地震时刻对应的应力值
      stress_rate_at_t : 地震时刻对应的应力变化率
    """
    # 计算应力变化率（利用向量化操作避免逐个循环）
    stress_rate = np.diff(Stress_cut) / np.diff(ttide_cut)

    # 找到离 t_ref 最近的索引
    i = np.argmin(np.abs(ttide_cut - t_ref))
    stress_at_t = Stress_cut[i]
    stress_rate_at_t = stress_rate[i] if i < len(stress_rate) else np.nan

    # 初始化返回变量
    ephase = np.nan
    maxtime = np.nan
    mintime1 = np.nan
    mintime2 = np.nan

    # 为相位插值构建 m 数组（单位为度，步长 0.5°）
    m_arr = np.arange(0, 180.5, 0.5)

    if stress_rate[i] >= 0:
        # ----- 正应力变化率情况 -----
        # 1. 找到从 i 开始，第一个使 stress_rate 变为负的索引 j
        candidates = np.where(stress_rate[i:] < 0)[0]
        j = i + candidates[0] if candidates.size > 0 else len(stress_rate)
        maxtime = ttide_cut[j - 1] if (j - 1) >= 0 else ttide_cut[i]

        # 2. 从 j 开始，找到第一个使 stress_rate 重新非负的索引 k
        candidates = np.where(stress_rate[j:] >= 0)[0]
        k = j + candidates[0] if candidates.size > 0 else len(stress_rate)
        mintime2 = ttide_cut[k - 1] if (k - 1) >= 0 else ttide_cut[i]

        # 3. 从起点到 i 内，寻找最后一个使 stress_rate 为负的索引 l
        candidates = np.where(stress_rate[:i+1] < 0)[0]
        if candidates.size > 0:
            l = candidates[-1]
            mintime1 = ttide_cut[l + 1] if (l + 1) < len(ttide_cut) else ttide_cut[i]
        else:
            mintime1 = ttide_cut[i]

        # 4. 向量化计算相位插值：求出每个 m 对应的 ephasetime
        ephasetime_arr = (maxtime - mintime1) / 180.0 * m_arr + mintime1
        valid = np.where(ttide_cut[i] <= ephasetime_arr)[0]
        if valid.size > 0:
            m_val = m_arr[valid[0]]
            ephase = -180 + m_val

    else:
        # ----- 负应力变化率情况 -----
        # 1. 从 i 开始，找到第一个使 stress_rate 变为非负的索引 j
        candidates = np.where(stress_rate[i:] >= 0)[0]
        j = i + candidates[0] if candidates.size > 0 else len(stress_rate)
        mintime2 = ttide_cut[j - 1] if (j - 1) >= 0 else ttide_cut[i]

        # 2. 从起点到 i 内，寻找最后一个使 stress_rate 为非负的索引 k
        candidates = np.where(stress_rate[:i+1] >= 0)[0]
        if candidates.size > 0:
            k = candidates[-1]
            maxtime = ttide_cut[k + 1] if (k + 1) < len(ttide_cut) else ttide_cut[i]
        else:
            maxtime = ttide_cut[i]

        # 3. 从 0 到 k 内，寻找最后一个使 stress_rate 为负的索引 l
        candidates = np.where(stress_rate[:k+1] < 0)[0]
        if candidates.size > 0:
            l = candidates[-1]
            mintime1 = ttide_cut[l + 1] if (l + 1) < len(ttide_cut) else ttide_cut[i]
        else:
            mintime1 = ttide_cut[i]

        # 4. 向量化计算相位插值
        ephasetime_arr = (mintime2 - maxtime) / 180.0 * m_arr + maxtime
        valid = np.where(ttide_cut[i] <= ephasetime_arr)[0]
        if valid.size > 0:
            m_val = m_arr[valid[0]]
            ephase = m_val

    # 转换相位单位（度 -> 弧度）
    if not np.isnan(ephase):
        ephase = np.deg2rad(ephase)
    else:
        ephase = np.nan

    return ephase, maxtime, mintime1, mintime2, stress_at_t, stress_rate_at_t


def fit_periodic_function(phi_degs, prob):
    """
    用最小二乘法拟合 R(φ) = a*cos(φ - φ0)
    其中 φ, φ0 均是度数，最终得到:
      - amplitude = a
      - phase_shift = φ0 (度)
      - model_func(φ) 可返回拟合后的 R(φ)
    参数:
      phi_degs : ndarray, 角度 (度)，形如 [0, 10, 20, ...]
      prob     : ndarray, 对应角度的观测值 (例如概率密度)
    返回:
      amplitude, phase_shift, model_func
    """
    # 1. 设计矩阵 G = [cos(φ), sin(φ)]
    G = np.column_stack((
        np.cos(np.deg2rad(phi_degs)),
        np.sin(np.deg2rad(phi_degs))
    ))

    # 2. 最小二乘求解 [A, B]
    #   R(φ) ~ A*cos(φ) + B*sin(φ)
    p, residuals, rank, s = np.linalg.lstsq(G, prob, rcond=None)
    A, B = p

    # 3. 转为复数 c = A - iB，便于同时求幅度和相位
    c = A - B*1j
    amplitude = np.abs(c)
    # 与 MATLAB 里: ph_shift = -angle(c, deg=True) 对齐
    phase_shift = -np.angle(c, deg=True)

    # 4. 构造预测函数
    def model_func(phi):
        # 输入 phi (度)，输出拟合的 R(φ)
        phi_rad = np.deg2rad(phi)
        return A*np.cos(phi_rad) + B*np.sin(phi_rad)

    return amplitude, phase_shift, model_func

#
#
# def modulation_phase(data, bin, PhStress, plot_option=None):
#     """
#     对应 MATLAB 的 ModulationPhase 函数
#
#     输入:
#       data     : 相位数据（单位：度）
#       bin      : 直方图分箱的中心值（与 MATLAB 中 hist 的第二个参数一致，要求均匀分布）
#       PhStress : 参考相位数据（单位：度）
#       gy       : 绘图选项，例如 'yes2' 或 'yes1'（可选）
#
#     输出:
#       PM_ra   : 归一化调制幅值
#       ph_shift: 相位偏移（度）
#       Po      : 均匀分布概率密度（每箱）
#       prob_1  : 原始直方图缩放后的值
#     """
#     data = np.asarray(data)
#     bin = np.asarray(bin)
#     PhStress = np.asarray(PhStress)
#
#     # MATLAB: Nnum = length(data);
#     Nnum = len(data)  # 虽然 MATLAB 中未在后续用到，但这里也计算一下
#
#     # MATLAB: [cout,ph1] = hist(data,bin);
#     # 为了与 MATLAB 的 hist 保持一致，我们需根据给定的 bin（箱中心）构造箱边界：
#     if len(bin) < 2:
#         raise ValueError("bin 数组至少需要包含两个元素以计算箱宽。")
#     bin_edges = np.empty(len(bin) + 1)
#     bin_edges[1:-1] = (bin[:-1] + bin[1:]) / 2
#     bin_edges[0] = bin[0] - (bin[1] - bin[0]) / 2
#     bin_edges[-1] = bin[-1] + (bin[-1] - bin[-2]) / 2
#     counts, _ = np.histogram(data, bins=bin_edges)
#     ph1 = bin.copy()  # 与 MATLAB 输出一致，ph1 就是 bin
#
#     # MATLAB: wbin = unique(diff(bin));
#     diffs = np.diff(bin)
#     unique_diffs = np.unique(diffs)
#     if unique_diffs.size == 1:
#         wbin = unique_diffs[0]
#     else:
#         raise ValueError("分箱中心的间距不均匀，不支持。")
#
#     # MATLAB: Prob_o = cout/length(data)/wbin;
#     Prob_o = counts / (len(data) * wbin)
#
#     # MATLAB: [cout2,~] = hist(PhStress,bin);
#     counts2, _ = np.histogram(PhStress, bins=bin_edges)
#
#     # MATLAB: Po = cout2/length(PhStress)/wbin;  % probability density
#     Po = counts2 / (len(PhStress) * wbin)
#
#     # MATLAB: G = [cosd(ph1)' sind(ph1)'];
#     G = np.column_stack((np.cos(np.deg2rad(ph1)), np.sin(np.deg2rad(ph1))))
#
#     # MATLAB: Prob = Prob_o./Po - 1;
#     Prob = Prob_o / Po - 1
#
#     # MATLAB: p = G(2:end-1,:)\Prob(2:end-1)';
#     if len(ph1) > 2:
#         G_sub = G[1:-1, :]  # MATLAB 中 2:end-1 对应 Python 的 1:-1
#         Prob_sub = Prob[1:-1]
#         p, _, _, _ = np.linalg.lstsq(G_sub, Prob_sub, rcond=None)
#     else:
#         p = np.zeros(2)
#
#     # MATLAB: c = p(1) - p(2)*sqrt(-1);
#     c = p[0] - p[1] * 1j
#
#     # MATLAB: Pm = abs(c);
#     Pm = np.abs(c)
#
#     # MATLAB: ph_shift = angle(c)*180/pi; ph_shift = -1*ph_shift;
#     ph_shift = -np.angle(c, deg=True)
#
#     # MATLAB: PM_ra = Pm/1;
#     PM_ra = Pm
#
#     # MATLAB: prob_1 = cout/(length(data)/length(cout));
#     prob_1 = counts / (len(data) / len(counts))
#
#     print(Prob,bin_edges)
#     phi_degs = (bin_edges[:-1] + bin_edges[1:]) / 2
#     a_fit, phi0_fit, model = fit_periodic_function(phi_degs ,Prob)
#     print(f"Fitted amplitude = {a_fit:.3f}")
#     print(f"Fitted phase_shift = {phi0_fit:.3f} deg")
#     print(f"Fitted amplitude = {PM_ra:.3f}")
#     print(f"Fitted phase_shift = {ph_shift:.3f} deg")
#     phi_dense = np.linspace(-180, 180, 30)
#
#
#     if plot_option is not None and plot_option.lower() == 'yes1':
#         vis.plot_phase_modulation(ph1, Prob_o, Po, G, p, wbin, PM_ra, ph_shift)
#     if plot_option is not None and plot_option.lower() == 'yes2':
#         vis.plot_phase_modulation_syn(ph1, Prob_o, Po, G, p, wbin, PM_ra, ph_shift,model,phi_dense)
#
#     return PM_ra, ph_shift, Po, prob_1



def modulation_phase(data, bin, PhStress, plot_option=None):
    """
    对应 MATLAB 的 ModulationPhase 函数

    输入:
      data     : 相位数据（单位：度）
      bin      : 直方图分箱的中心值（与 MATLAB 中 hist 的第二个参数一致，要求均匀分布）
      PhStress : 参考相位数据（单位：度）
      gy       : 绘图选项，例如 'yes2' 或 'yes1'（可选）

    输出:
      PM_ra   : 归一化调制幅值
      ph_shift: 相位偏移（度）
      Po      : 均匀分布概率密度（每箱）
      prob_1  : 原始直方图缩放后的值
    """
    data = np.asarray(data)
    bin = np.asarray(bin)
    PhStress = np.asarray(PhStress)


    # MATLAB: [cout,ph1] = hist(data,bin);
    # 为了与 MATLAB 的 hist 保持一致，我们需根据给定的 bin（箱中心）构造箱边界：
    if len(bin) < 2:
        raise ValueError("bin 数组至少需要包含两个元素以计算箱宽。")
    bin_edges = np.empty(len(bin) + 1)
    bin_edges[1:-1] = (bin[:-1] + bin[1:]) / 2
    bin_edges[0] = bin[0] - (bin[1] - bin[0]) / 2
    bin_edges[-1] = bin[-1] + (bin[-1] - bin[-2]) / 2
    counts, cc1 = np.histogram(data, bins=bin_edges)
    ph1 = bin.copy()  # 与 MATLAB 输出一致，ph1 就是 bin

    # MATLAB: wbin = unique(diff(bin));
    diffs = np.diff(bin)
    unique_diffs = np.unique(diffs)
    if unique_diffs.size == 1:
        wbin = unique_diffs[0]
    else:
        raise ValueError("分箱中心的间距不均匀，不支持。")

    # MATLAB: Prob_o = cout/length(data)/wbin;
    Prob_o = counts / (len(data) * wbin)

    # MATLAB: [cout2,~] = hist(PhStress,bin);
    counts2, _ = np.histogram(PhStress, bins=bin_edges)

    # MATLAB: Po = cout2/length(PhStress)/wbin;  % probability density
    Po = counts2 / (len(PhStress) * wbin)

    # MATLAB: G = [cosd(ph1)' sind(ph1)'];
    G = np.column_stack((np.cos(np.deg2rad(ph1)), np.sin(np.deg2rad(ph1))))

    # MATLAB: Prob = Prob_o./Po - 1;
    Prob = Prob_o / Po - 1

    phi_degs = (bin_edges[:-1] + bin_edges[1:]) / 2
    # print(Prob, phi_degs )
    # print(len(bin_edges),len(phi_degs))
    a_fit, phi0_fit, model = fit_periodic_function(phi_degs ,Prob)
    print(f"Fitted amplitude = {a_fit:.3f}")
    print(f"Fitted phase_shift = {phi0_fit:.3f} deg")
    phi_dense = np.linspace(-180, 180, 30)
    ph_shift  = phi0_fit
    PM_ra = a_fit

    if plot_option is not None and plot_option.lower() == 'yes1':
        # vis.plot_phase_modulation(ph1, Prob_o, Po, G, p, wbin, PM_ra, ph_shift)
        vis.plot_phase_modulation(ph1, Prob_o, Po, wbin, PM_ra, ph_shift,model,phi_dense)

    return PM_ra, ph_shift, Po, Prob




def modulation_stress(Stress_AM_bk, Stress_AM, initial_guess, plot_option):
    """
    对应 MATLAB 的 ModulationStress 函数
    Stress_AM_S_bk: 参考应力数据（一维数组）
    Stress_AM_S: 观测应力数据（一维数组）
    initial_guess: 优化初始参数，如 [a0, C0]
    plot_option: 'yes' 绘图，'no' 不绘图
    返回： a_estimated, C_estimated, delta_a, delta_c
    """

    # 1. 分箱
    edges = np.linspace(np.min(Stress_AM), np.max(Stress_AM), 20)
    counts, _ = np.histogram(Stress_AM, bins=edges)
    counts_bk, _ = np.histogram(Stress_AM_bk, bins=edges)
    event_rate = counts / counts_bk

    # 2. 过滤
    nonzero_idx = event_rate>0
    # counts = counts[nonzero_idx]
    # counts_bk = counts_bk[nonzero_idx]
    shear_stress = (edges[:-1] + edges[1:]) / 2
    shear_stress = shear_stress[nonzero_idx]

    # 3. 转换为 kPa
    shear_stress_kPa = shear_stress / 1000.0

    counts = counts[nonzero_idx]
    counts_bk = counts_bk[nonzero_idx]
    event_rate = counts / counts_bk
    # print(event_rate)
    # print(counts)
    # print(counts_bk)

    #
    # def neg_log_likelihood(params):
    #     a, C = params
    #     model_rate = C * np.exp(a * shear_stress_kPa)
    #     model_rate = np.where(model_rate <= 0, 1e-25, model_rate)
    #     return -np.sum(counts * np.log(model_rate) - counts_bk * model_rate)
    #
    # options = {
    #     # 对应 MATLAB 的 'Display'='off'
    #     'disp': False,
    #     # 对应 MATLAB 的 'MaxIterations' = 400
    #     'maxiter': 400,
    #     # 对应 MATLAB 的 'OptimalityTolerance' = 1e-6
    #     # SciPy 中 BFGS 的梯度收敛阈值
    #     'gtol': 1e-6,
    #     # 对应 MATLAB 的 'FiniteDifferenceStepSize' = sqrt(eps)
    #     # 在数值近似梯度时使用的步长
    #     'eps': np.sqrt(np.finfo(float).eps),
    #     # 如果想让函数值或步长变化小于某阈值就停止，可以再加一个全局收敛判定
    #     'tol': 1e-6,   # 这会作为一个总的收敛容忍度
    # }
    #
    # res = minimize(neg_log_likelihood, initial_guess, method='BFGS',options=options)
    # print(res)
    # a_estimated, C_estimated = res.x
    # hess_inv = res.hess_inv
    #
    # # 计算误差
    # std_error_a = np.sqrt(hess_inv[0, 0])
    # std_error_C = np.sqrt(hess_inv[1, 1])
    # delta_a = 2 * std_error_a
    # delta_c = 2 * std_error_C

    def neg_log_likelihood(params):
        a, C = params
        # 为避免 log(0) 或负值，通常加个下限
        model_rate = C * np.exp(a * shear_stress_kPa)
        model_rate = np.where(model_rate <= 0, 1e-60, model_rate)
        return -np.sum(counts * np.log(model_rate) - counts_bk * model_rate)

    def grad_neg_log_likelihood(params):
        a, C = params
        # 同步计算 model_rate
        model_rate = C * np.exp(a * shear_stress_kPa)
        model_rate = np.where(model_rate <= 0, 1e-60, model_rate)
        # d(nll)/da
        grad_a = np.sum(
            shear_stress_kPa * (counts_bk * model_rate - counts)
        )
        # d(nll)/dC
        grad_C = np.sum(
            counts_bk * np.exp(a * shear_stress_kPa) - counts / C
        )
        return np.array([grad_a, grad_C])

    def hess_neg_log_likelihood(params):
        a, C = params
        # 计算二阶导数
        # 先做一些中间量，便于重复利用
        exp_term = np.exp(a * shear_stress_kPa)  # e^{a x_i}
        model_rate = C * exp_term  # C e^{a x_i}

        # H[0,0] = d^2(nll)/(da^2)
        h_aa = np.sum(shear_stress_kPa ** 2 * counts_bk * model_rate)

        # H[1,1] = d^2(nll)/(dC^2)
        h_CC = np.sum(counts / (C ** 2))

        # H[0,1] = H[1,0] = d^2(nll)/(da dC)
        h_aC = np.sum(shear_stress_kPa * counts_bk * exp_term)

        # 拼成 2x2 矩阵
        return np.array([
            [h_aa, h_aC],
            [h_aC, h_CC]
        ])


    options = {
        'disp': False,  # 是否打印优化信息
        'maxiter': 400,  # 最大迭代次数
        'xtol': 1e-6,  # 或者 'tol'、'gtol' 等可做精度调节
    }

    res = minimize(
        fun=neg_log_likelihood,
        x0=initial_guess,
        method='Newton-CG',
        jac=grad_neg_log_likelihood,
        hess=hess_neg_log_likelihood,
        options=options
    )
    # print(res)
    a_estimated, C_estimated = res.x

    H = hess_neg_log_likelihood(res.x)
    cov_matrix = np.linalg.inv(H)

    # print(cov_matrix)
    # print(x)

    # 提取对角线元素的方差并开方
    std_error_a = np.sqrt(cov_matrix[0, 0])
    std_error_C = np.sqrt(cov_matrix[1, 1])
    # 计算 95% 置信区间（常用近似 2 倍标准误差）
    delta_a = 2 * std_error_a
    delta_c = 2 * std_error_C


    # print("Python hess_inv:\n", hess_inv)
    print("Estimated a:", a_estimated, "with delta:", delta_a)
    print("Estimated C:", C_estimated, "with delta:", delta_c)
    if plot_option == 'yes':
        # 调用单独的绘图函数
        fig, ax1, ax2 = vis.plot_modulation_stress(
            Stress_AM_bk, Stress_AM, shear_stress, shear_stress_kPa,
            event_rate, a_estimated, C_estimated, delta_a, delta_c
        )

        # 注意：这里不调用 plt.show()，由外部统一管理

    return a_estimated, C_estimated, delta_a, delta_c, Stress_AM_bk, Stress_AM, shear_stress, shear_stress_kPa, event_rate



def schuster_test(pha):
    """Schuster 检验的 p-value 计算函数
    参数:
    pha : ndarray
        以弧度为单位的相位数组

    返回:
    p_value : float
        Schuster 检验的 p-value
    """
    # 转换为角度
    pha_deg = np.degrees(pha)

    # 计算 R_x 和 R_y
    R_x = np.sum(np.cos(np.radians(pha_deg)))
    R_y = np.sum(np.sin(np.radians(pha_deg)))

    # 计算 R 和 Z
    R = np.sqrt(R_x ** 2 + R_y ** 2)
    N = len(pha)
    Z = (R ** 2) / N

    # 计算 p-value
    p_value = np.exp(-Z)
    return p_value



def grid_modulation_results_2D(Allpha_obs, Allpha_ref, bin_phase, initial_guess, grid_spacing, threshold):
    """
    Performs 2D gridding of the data where each grid cell has a fixed size of grid_spacing degrees.

    Parameters:
        Allpha_obs (ndarray): Observation data with columns representing [?, phase, stress, ..., lat, lon, ...].
                             Assumes that:
                             - Column 2 (index 1) is phase.
                             - Column 3 (index 2) is stress.
                             - Column 11 (index 10) is latitude.
                             - Column 12 (index 11) is longitude.
        Allpha_ref (ndarray): Reference data with similar formatting. Uses:
                             - Column 2 (index 1) for phase.
                             - Column 3 (index 2) for stress.
        bin_phase: Parameter for the ModulationPhase function.
        initial_guess: Initial guess for the ModulationStress function.
        grid_spacing (float): Grid size in degrees (e.g., 0.2).
        threshold (int): Minimum number of samples required in a grid cell to perform calculations.

    Returns:
        grid_centers (ndarray): An array of shape (n, 7) where each row contains:
            [lon_center, lat_center, a_val, p_val, pm_val, pm_pha, number_eq]
    """
    # Extract longitude and latitude from observation data
    # Note: Adjust indices to match your data structure (MATLAB columns 12 and 11 correspond to Python indices 11 and 10)
    lon = Allpha_obs[:, 11]
    lat = Allpha_obs[:, 10]

    # Define grid edges in longitude and latitude directions
    lon_edges = np.arange(np.floor(np.min(lon)), np.ceil(np.max(lon)) + grid_spacing, grid_spacing)
    lat_edges = np.arange(np.floor(np.min(lat)), np.ceil(np.max(lat)) + grid_spacing, grid_spacing)

    # Calculate grid centers
    lon_centers = (lon_edges[:-1] + lon_edges[1:]) / 2
    lat_centers = (lat_edges[:-1] + lat_edges[1:]) / 2

    # Determine the number of grid cells in each direction
    num_bins_lon = len(lon_centers)
    num_bins_lat = len(lat_centers)

    # List to collect grid center results; each row: [lon, lat, a, p, pm, pm_pha, count]
    grid_centers_vol = []

    # Loop over grid cells (lon index first, then lat index)
    for i in range(num_bins_lon):
        for j in range(num_bins_lat):
            # Identify the samples within the current grid cell
            idx = ((lon >= lon_edges[i]) & (lon < lon_edges[i + 1]) &
                   (lat >= lat_edges[j]) & (lat < lat_edges[j + 1]))
            count = np.sum(idx)

            # Skip grid cells with insufficient samples
            if count < threshold:
                continue

            # Extract phase and stress data from the observation data.
            # In MATLAB, column 2 is phase and column 3 is stress.
            phase_obs = Allpha_obs[idx, 1]
            stress_obs = Allpha_obs[idx, 2]

            # Calculate modulation parameters.
            # Convert phases from radians to degrees.
            # Note: The string 'no' is passed as an option (presumably to suppress plotting or messages).

            pm_val, pm_pha, _, _ = modulation_phase(phase_obs * 180 / np.pi,bin_phase,Allpha_ref[:, 1] * 180 / np.pi, 'no')

            # Calculate stress modulation parameter 'a'
            a_val,_, _,_,_, _,_ ,_, _ = modulation_stress(Allpha_ref[:, 2], stress_obs, initial_guess, 'no')

            # Calculate p-value using the Schuster Test
            p_val = schuster_test(phase_obs)

            # Determine the grid cell center coordinates
            lon_center = lon_centers[i]
            lat_center = lat_centers[j]

            # Append the results to the list
            grid_centers_vol.append([lon_center, lat_center, a_val, p_val, pm_val, pm_pha, count])

    # Convert the result to a NumPy array and return
    grid_centers = np.array(grid_centers_vol)
    return grid_centers
