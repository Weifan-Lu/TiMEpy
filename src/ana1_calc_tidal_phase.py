import numpy as np
import datetime
import matplotlib.dates as mdates
# 假设 date2num 已经从 matplotlib.dates 导入，例如：
from matplotlib.dates import date2num
import main_function as mf
import visualization as vis

def ana1_calc_tidal_phase(config, opt):
    """
    根据参数字典 params 和模式 opt 计算潮汐相位，并将结果保存为 TXT 文件。

    params 字典中需包含（部分参数均以 MATLAB datenum 格式给出，单位：天）：
      - data_select_decluster : 地震目录文件（若 opt=='1' 使用该文件）
      - data_select           : 地震目录文件（若 opt!='1' 使用该文件）
      - stress_to_mat         : 包含应力数据的 .mat 文件路径（其中变量 normal_stress, shear_stress, volumetric_strain）
      - t_stress_start        : 应力数据起始时间（字符串格式 'YYYY-MM-DD' 或 datenum 数值）
      - t_sample              : 采样间隔（秒）
      - miu                   : 摩擦系数
      - start_time            : 分析起始时间（datenum 数值）
      - end_time              : 分析结束时间（datenum 数值）
      - t_search              : 搜索窗口（单位：天，注意：MATLAB 代码中用 t_search*86400/t_sample 计算样本数）
      - ref_interval          : 均匀参考时间间隔（天）
      - phase_stress_obs      : 保存观测结果的 TXT 文件名前缀（会自动附加 _N.txt, _S.txt, 等）
      - phase_stress_ref      : 保存参考结果的 TXT 文件名前缀
    opt 为 '1' 时使用 data_select_decluster，否则使用 data_select。
    """
    # --------------------
    # 1. 读取地震目录数据
    if opt == '1':
        data = np.loadtxt(config.data_select_decluster)
    else:
        data = np.loadtxt(config.data_select)
    # 假设 data 每行的列依次为：
    # year, month, day, hour, minute, sec, (其他)，lat, lon, depth
    years = data[:, 0].astype(int)
    months = data[:, 1].astype(int)
    days = data[:, 2].astype(int)
    hours = data[:, 3].astype(int)
    minutes = data[:, 4].astype(int)
    secs = data[:, 5]
    # 构造 datetime 对象列表
    eq_datetimes = [datetime.datetime(y, m, d, h, mi, int(s))
                    for y, m, d, h, mi, s in zip(years, months, days, hours, minutes, secs)]
    # 转换为 datenum（matplotlib 的 date2num 与 MATLAB datenum 格式基本兼容）
    t_decluster = mdates.date2num(eq_datetimes)
    print(t_decluster)

    # --------------------
    # 2. 读取应力数据
    # 只读取需要的列（假设文件中无标题行，如果有标题行则可增加 skiprows=1）
    data_stress = np.loadtxt(config.output_stress_txt, usecols=(3, 4, 5))

    # 根据写入格式，索引对应：
    volumetric_strain = data_stress[:, 0]
    shear_stress = data_stress[:, 1]
    normal_stress = data_stress[:, 2]
    miu = config.miu
    StressN = normal_stress
    StressS = shear_stress
    StressVol = volumetric_strain
    StressCFS = shear_stress + miu * normal_stress


    # --------------------
    # 3. 构造应力数据的时间序列 ttide
    # t_stress_start 可为字符串或数值
    if isinstance(config.t_stress_start, str):
        t0_dt = datetime.datetime.strptime(config.t_stress_start, '%Y-%m-%d')
        t_stress_start = mdates.date2num(t0_dt)
    else:
        t_stress_start = config.t_stress_start
    t_sample = config.t_sample  # 单位：秒
    n_stress = len(normal_stress)
    ttide = t_stress_start + np.arange(n_stress) * t_sample / 86400.0  # 转换为天
    t_tide = ttide.copy()


    # --------------------
    # 4. 对目录数据（obs）计算潮汐相位
    t_start = config.start_time  # datenum 数值
    t_end = config.end_time
    print('start time:',t_start,'  end time: ',t_end,' len_stress: ',len(StressN),'len_time: ',len(ttide))
    # 如果 t_start 和 t_end 是 datetime 对象，则转换为 datenum 格式（float）
    if isinstance(t_start, datetime.datetime):
        t_start = mdates.date2num(t_start)
    if isinstance(t_end, datetime.datetime):
        t_end = mdates.date2num(t_end)
    indices = np.where((t_decluster >= t_start) & (t_decluster <= t_end))[0]

    obs_list_N = []
    obs_list_S = []
    obs_list_CFS = []
    obs_list_Vol = []

    incat = int(round(config.t_search * 86400 / t_sample))
    print(indices)

    for ix in indices:
        t_cut = t_decluster[ix]
        Index = np.argmin(np.abs(t_tide - t_cut))
        range_indices = np.arange(Index - incat, Index + incat + 1)
        # print(range_indices)
        # 检查索引范围是否越界
        if range_indices[0] < 0 or range_indices[-1] >= len(StressN):
            # print('wrong')
            continue
        StressN_cut = StressN[range_indices]
        StressS_cut = StressS[range_indices]
        StressVol_cut = StressVol[range_indices]
        StressCFS_cut = StressCFS[range_indices]
        ttide_cut = t_tide[range_indices]

        # 分别计算潮汐相位等信息（这里只用其中的相位、应力值和应力变化率）
        phaN, _, _, _, stress_at_t_N, stress_rate_at_t_N = mf.calculate_tidal_phase(ttide_cut, StressN_cut, t_cut)
        phaS, _, _, _, stress_at_t_S, stress_rate_at_t_S = mf.calculate_tidal_phase(ttide_cut, StressS_cut, t_cut)
        phaCFS, _, _, _, stress_at_t_CFS, stress_rate_at_t_CFS = mf.calculate_tidal_phase(ttide_cut, StressCFS_cut, t_cut)
        phaVol, _, _, _, stress_at_t_Vol, stress_rate_at_t_Vol = mf.calculate_tidal_phase(ttide_cut, StressVol_cut, t_cut)
        dt = mdates.num2date(t_cut)
        obs_list_N.append([t_cut, phaN, stress_at_t_N, stress_rate_at_t_N,dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second])
        obs_list_S.append([t_cut, phaS, stress_at_t_S, stress_rate_at_t_S])
        obs_list_CFS.append([t_cut, phaCFS, stress_at_t_CFS, stress_rate_at_t_CFS])
        obs_list_Vol.append([t_cut, phaVol, stress_at_t_Vol, stress_rate_at_t_Vol])
    AllphaN = np.array(obs_list_N)
    AllphaS = np.array(obs_list_S)
    AllphaCFS = np.array(obs_list_CFS)
    AllphaVol = np.array(obs_list_Vol)
    print('len obs phas:',len(AllphaVol))

    vis.plot_stress_earthquake(AllphaVol[:, 0],AllphaVol[:, 2], AllphaN[:, 2], AllphaS[:, 2],AllphaCFS[:, 2],t_tide,StressVol,StressN,StressS,StressCFS,config)


    # --------------------
    # 5. 生成均匀参考时间序列并计算潮汐相位（ref）
    ref_interval = config.ref_interval  # 单位：天
    n_samples = int(round((t_end - t_start) / ref_interval)) + 1
    t_ref = np.linspace(t_start, t_end, n_samples)

    ref_list_N = []
    ref_list_S = []
    ref_list_CFS = []
    ref_list_Vol = []
    for t_cut in t_ref:
        Index = np.argmin(np.abs(t_tide - t_cut))
        range_indices = np.arange(Index - incat, Index + incat + 1)
        if range_indices[0] < 0 or range_indices[-1] >= len(StressN):
            continue
        StressN_cut = StressN[range_indices]
        StressS_cut = StressS[range_indices]
        StressVol_cut = StressVol[range_indices]
        StressCFS_cut = StressCFS[range_indices]
        ttide_cut = t_tide[range_indices]

        phaN, _, _, _, stress_at_t_N, stress_rate_at_t_N = mf.calculate_tidal_phase(ttide_cut, StressN_cut, t_cut)
        phaS, _, _, _, stress_at_t_S, stress_rate_at_t_S = mf.calculate_tidal_phase(ttide_cut, StressS_cut, t_cut)
        phaCFS, _, _, _, stress_at_t_CFS, stress_rate_at_t_CFS = mf.calculate_tidal_phase(ttide_cut, StressCFS_cut, t_cut)
        phaVol, _, _, _, stress_at_t_Vol, stress_rate_at_t_Vol = mf.calculate_tidal_phase(ttide_cut, StressVol_cut, t_cut)

        ref_list_N.append([t_cut, phaN, stress_at_t_N, stress_rate_at_t_N])
        ref_list_S.append([t_cut, phaS, stress_at_t_S, stress_rate_at_t_S])
        ref_list_CFS.append([t_cut, phaCFS, stress_at_t_CFS, stress_rate_at_t_CFS])
        ref_list_Vol.append([t_cut, phaVol, stress_at_t_Vol, stress_rate_at_t_Vol])
    AllphaN_ref = np.array(ref_list_N)
    AllphaS_ref = np.array(ref_list_S)
    AllphaCFS_ref = np.array(ref_list_CFS)
    AllphaVol_ref = np.array(ref_list_Vol)

    # --------------------
    # 6. 保存结果为 TXT 文件
    header = "t (datenum)\tphase (rad)\tstress_at_t\tstress_rate_at_t"
    np.savetxt(config.phase_stress_obs + '_N.txt', AllphaN, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_obs + '_S.txt', AllphaS, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_obs + '_CFS.txt', AllphaCFS, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_obs + '_Vol.txt', AllphaVol, fmt='%.8f', delimiter='\t', header=header)

    np.savetxt(config.phase_stress_ref + '_N.txt', AllphaN_ref, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_ref + '_S.txt', AllphaS_ref, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_ref + '_CFS.txt', AllphaCFS_ref, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_ref + '_Vol.txt', AllphaVol_ref, fmt='%.8f', delimiter='\t', header=header)

    print("潮汐相位结果已保存为 TXT 文件。")