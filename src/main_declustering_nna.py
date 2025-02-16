import numpy as np
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from scipy.ndimage import uniform_filter


def haversine_distance(lat1, lon1, lat2, lon2, radius=6378.1):
    """
    计算大圆距离（单位：公里），支持标量或数组输入，角度单位为度。
    """
    # 转换为弧度
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    # 球面余弦公式
    cos_val = np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2)
    cos_val = np.clip(cos_val, -1, 1)
    d = np.arccos(cos_val)
    return d * radius


def smooth2a(A, Nr, Nc):
    """
    对二维数组 A 进行移动平均平滑，窗口大小为 (Nr, Nc)
    """
    return uniform_filter(A, size=(Nr, Nc))


def decluster_nna(start_time,data,df_val,b,p,q,eta0):
    """
    根据控制参数处理地震目录数据，封装了计算 Nij, Tij, Rij，构造 BgLink/StLink、家庭识别等流程。

    参数：
        params (dict): 控制参数字典，必须包含如下键：
            - 'start_time': （可选）开始时间（此示例中未使用）
            - 'data_select': 数据文件路径（当未传入 data 时使用）
            - 'df': 距离幂指数（例如 1.6）
            - 'b': b 值参数
            - 'p': p 值参数
            - 'q': q 值参数
            - 'eta0': 用于分类的阈值（对 LinkDis 第三列取对数后比较）
        data (np.ndarray, optional): 数据数组，若传入则忽略 data_select。数据要求形状为 (N, 10)，各列依次为：
            [year, month, day, hour, minute, second, lat, lon, dep, mag]

    返回：
        dict: 包含处理结果的字典，键包括：
            - 'Declustered_Catalog': 去相关目录，每行数据为原 data 对应行后附加 t/365（年）
            - 'Dataobj': 后处理后的 Dataobj 数组（第一、二列取了对数，第三列为 min_nj）
            - 'LinkDis': LinkDis 数组（第三列已取对数）
            - 'BgLink': BgLink 数组
            - 'StLink': StLink 数组
            - 'family_weak': weak family 数组
            - 'Mainshock': 主震家庭字典，每个键对应一个家庭的事件索引数组
            - 'All_Mainshock': 所有主震事件的索引列表
            - 'All_Foreshock': 前震组列表（每组为二维数组）
            - 'All_Aftershock': 后震组列表（每组为二维数组）
            - 'All_Family': 合并后的家庭事件索引数组
            - 'All_Single': 单一事件（孤立事件）的集合
            - 'z': Dataobj 的第三列（min_nj）
    """

    # 提取各列数据（注意：Python 索引从 0 开始）
    year = data[:, 0].astype(int)
    month = data[:, 1].astype(int)
    day = data[:, 2].astype(int)
    hour = data[:, 3].astype(int)
    minute = data[:, 4].astype(int)
    sec = data[:, 5]  # 秒可能包含小数部分
    lat = data[:, 6]
    lon = data[:, 7]
    dep = data[:, 8]
    mag = data[:, 9]

    # 将时间信息转换为 datetime 对象，再转换为 matplotlib 日期数（单位：天）
    t_list = [datetime(y, m, d, h, mi, int(s))
              for y, m, d, h, mi, s in zip(year, month, day, hour, minute, sec)]
    t = mdates.date2num(t_list)

    # 数据行数
    len_data = data.shape[0]

    # ======== 计算 Nij, Tij, Rij, Dataobj, LinkDis ========
    Dataobj = []  # 存储 [min_tj, min_rj, min_nj]
    LinkDis = []  # 存储 [当前事件索引, 最近邻事件索引, min_nj]
    for ii in range(1, len_data):  # ii 从 1 到 len_data-1
        dt = (t[ii] - t[:ii]) / 365.0  # 时间差（年）
        dr = haversine_distance(lat[ii], lon[ii], lat[:ii], lon[:ii])
        dm = mag[:ii]
        Nij = dt * (dr ** df_val) * (10 ** (-b * dm))
        Tij = dt * (10 ** (-q * b * dm))
        Rij = (dr ** df_val) * (10 ** (-p * b * dm))
        # 将值为 0 的元素替换为 nan
        Nij = np.where(Nij == 0, np.nan, Nij)
        Tij = np.where(Tij == 0, np.nan, Tij)
        Rij = np.where(Rij == 0, np.nan, Rij)
        if np.all(np.isnan(Nij)):
            continue
        # 找到 Nij 中的最小值及其索引（忽略 nan）
        xnumx = np.nanargmin(Nij)
        min_nj = Nij[xnumx]
        min_tj = Tij[xnumx]
        min_rj = Rij[xnumx]
        Dataobj.append([min_tj, min_rj, min_nj])
        LinkDis.append([ii, xnumx, min_nj])
    Dataobj = np.array(Dataobj)
    LinkDis = np.array(LinkDis)

    # 对 LinkDis 第三列取对数
    if LinkDis.size > 0:
        LinkDis[:, 2] = np.log10(LinkDis[:, 2])

    # ======== 构造 BgLink 与 StLink ============
    BgLink = []
    StLink = []
    for i in range(LinkDis.shape[0]):
        # 若 log10(min_nj) > eta0，则归为 BgLink，否则归为 StLink
        if LinkDis[i, 2] > eta0:
            BgLink.append([int(LinkDis[i, 0]), int(LinkDis[i, 1]), LinkDis[i, 2],
                           t[int(LinkDis[i, 0])], t[int(LinkDis[i, 1])],
                           lat[int(LinkDis[i, 0])], lat[int(LinkDis[i, 1])]])
        if LinkDis[i, 2] < eta0:
            StLink.append([int(LinkDis[i, 0]), int(LinkDis[i, 1]), LinkDis[i, 2],
                           t[int(LinkDis[i, 0])], t[int(LinkDis[i, 1])],
                           lat[int(LinkDis[i, 0])], lat[int(LinkDis[i, 1])]])
    BgLink = np.array(BgLink)
    StLink = np.array(StLink)

    # ======== 识别 family_weak ============
    family_weak = []
    if BgLink.size > 0 and StLink.size > 0:
        for row in BgLink:
            if int(row[0]) in StLink[:, 1].astype(int):
                family_weak.append([int(row[0]), 1])
    family_weak = np.array(family_weak) if len(family_weak) > 0 else np.empty((0, 2))

    # ======== 构造家庭索引 ============
    num_st = StLink.shape[0]
    Family_idx = np.arange(num_st)  # 初始家庭编号
    for i in range(1, num_st):
        for j in range(i):
            if (StLink[j, 0] == StLink[i, 1]) or (StLink[j, 1] == StLink[i, 1]):
                Family_idx[i] = Family_idx[j]
                break
    Fm = np.unique(Family_idx)

    # ======== 构造 Mainshock 家族 ============
    Mainshock = {}  # 使用字典存储，每个键对应一个家庭
    for i in range(num_st):
        fam = Family_idx[i]
        if fam not in Mainshock:
            Mainshock[fam] = []
        Mainshock[fam].append(StLink[i, 0:2])
    Mainshock = {k: np.array(v) for k, v in Mainshock.items() if len(v) > 0}

    All_Mainshock = []
    All_Foreshock = []
    All_Aftershock = []
    for key in Mainshock:
        aa = Mainshock[key]
        # 将家庭中所有事件索引合并后去重
        AA = np.unique(np.concatenate((aa[:, 0], aa[:, 1])).astype(int))
        # 选择震级最大的事件为主震
        mags = mag[AA]
        num_max = np.argmax(mags)
        mainshock_idx = AA[num_max]
        All_Mainshock.append(mainshock_idx)
        # 构造一个二维数组，第二列均为主震索引
        AB = np.column_stack((AA, np.full(AA.shape, mainshock_idx)))
        if num_max == len(AB) - 1:
            Aftershock = np.empty((0, 2), dtype=int)
            Foreshock = AB[:num_max, :]
        elif num_max == 0:
            Foreshock = np.empty((0, 2), dtype=int)
            Aftershock = AB[1:, :]
        else:
            Foreshock = AB[:num_max, :]
            Aftershock = AB[num_max + 1:, :]
        if Foreshock.size > 0:
            All_Foreshock.append(Foreshock)
        if Aftershock.size > 0:
            All_Aftershock.append(Aftershock)

    # 合并各家庭中的事件索引
    family_list = []
    if All_Foreshock:
        for group in All_Foreshock:
            family_list.append(group[:, 0])
    if All_Aftershock:
        for group in All_Aftershock:
            family_list.append(group[:, 0])
    if All_Mainshock:
        family_list.append(np.array(All_Mainshock))
    if family_list:
        All_Family = np.unique(np.concatenate(family_list))
    else:
        All_Family = np.array([])

    # ======== 单一事件 ============
    # 注意：MATLAB 中的 1 对应 Python 的 0（第一个事件）
    single_candidates = set(BgLink[:, 0].astype(int)) if BgLink.size > 0 else set()
    weak_family = set(family_weak[:, 0].astype(int)) if family_weak.size > 0 else set()
    All_Single = {0} | (single_candidates - weak_family)

    # ======== 去相关事件 ============
    All_decluster = np.unique(np.concatenate((np.array(list(All_Single)), np.array(All_Mainshock))))
    Declustered_Catalog = []
    for idx in All_decluster:
        Declustered_Catalog.append(np.concatenate((data[int(idx), :], [t[int(idx)] / 365.0])))
    Declustered_Catalog = np.array(Declustered_Catalog)

    # ======== Dataobj 后处理 ============
    if Dataobj.size > 0:
        Dataobj[:, 0] = np.log10(Dataobj[:, 0])
        Dataobj[:, 1] = np.log10(Dataobj[:, 1])
        z = Dataobj[:, 2]
    else:
        z = np.array([])

    return {
        'Declustered_Catalog': Declustered_Catalog,
        'Dataobj': Dataobj,
        'LinkDis': LinkDis,
        'BgLink': BgLink,
        'StLink': StLink,
        'family_weak': family_weak,
        'Mainshock': Mainshock,
        'All_Mainshock': All_Mainshock,
        'All_Foreshock': All_Foreshock,
        'All_Aftershock': All_Aftershock,
        'All_Family': All_Family,
        'All_Single': All_Single,
        'z': z,
        't': t
    }
    # fig3, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
    # ax1.plot(t, data[:, 6], '.', markersize=4)
    # ax1.set_ylabel('Lat')
    # ax1.xaxis_date()
    # fig3.autofmt_xdate()
    # ax1.tick_params(labelsize=16)
    #
    # # 根据去相关目录前 6 列重构时间
    # years_out = Declustered_Catalog[:, 0].astype(int)
    # months_out = Declustered_Catalog[:, 1].astype(int)
    # days_out = Declustered_Catalog[:, 2].astype(int)
    # hours_out = Declustered_Catalog[:, 3].astype(int)
    # minutes_out = Declustered_Catalog[:, 4].astype(int)
    # secs_out = Declustered_Catalog[:, 5]
    # t1_list = []
    # for i in range(len(years_out)):
    #     t1_list.append(datetime(years_out[i], months_out[i], days_out[i],
    #                             hours_out[i], minutes_out[i], int(secs_out[i])))
    # t1 = mdates.date2num(t1_list)
    # ax2.plot(t1, Declustered_Catalog[:, 6], '.', markersize=4)
    # ax2.set_ylabel('Lat')
    # ax2.xaxis_date()
    # fig3.autofmt_xdate()
    # ax2.tick_params(labelsize=16)
    # fig3.tight_layout()
    # fig3.savefig(config.decluster_data_fig_lat, dpi=300)
    #
    # # ----- Figure 5: 二维直方图及平滑图像 -----
    # bin_width = 0.1
    # num_bins = int(round(10 / bin_width))
    # X = Dataobj[:, 0]
    # Y = Dataobj[:, 1]
    # A = np.zeros((num_bins, num_bins))
    # # 计算每个网格中事件数
    # for i in range(num_bins):
    #     for j in range(num_bins):
    #         x_min = -9 - bin_width + bin_width * (i + 1)
    #         x_max = -9 + bin_width * (i + 1)
    #         y_min = -6 - bin_width + bin_width * (j + 1)
    #         y_max = -6 + bin_width * (j + 1)
    #         A[j, i] = np.sum((X >= x_min) & (X < x_max) & (Y > y_min) & (Y < y_max))
    # x_extent = [-9, -9 + bin_width * num_bins]
    # y_extent = [-6, -6 + bin_width * num_bins]
    #
    # fig5, ax = plt.subplots(figsize=(6, 4))
    # im = ax.imshow(A, extent=(x_extent[0], x_extent[1], y_extent[0], y_extent[1]),
    #                origin='lower', cmap='jet')
    # plt.colorbar(im, ax=ax)
    # # 平滑处理
    # Nr, Nc = 2, 2
    # matrixOut = mdnna.smooth2a(A, Nr, Nc)
    # ax.imshow(matrixOut, extent=(x_extent[0], x_extent[1], y_extent[0], y_extent[1]),
    #           origin='lower', cmap='jet')
    # ax.set_xlim([-7, -1])
    # ax.set_ylim([-5, 3])
    # # 绘制绿色直线：T + R = constant，即 R = eta0 - T
    # T_line = np.arange(-10, 3, 1)
    # R_line = eta0 - T_line
    # ax.plot(T_line, R_line, color='white', linewidth=1)
    # ax.set_xlabel('Rescaled time, T')
    # ax.set_ylabel('Rescaled distance, R')
    # ax.set_xticks([-7, -6, -5, -4, -3, -2, -1])
    # ax.set_yticks([-5, -4, -3, -2, -1, 0, 1, 2, 3])
    # ax.set_xticklabels(['10^{-7}', '10^{-6}', '10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}'])
    # ax.set_yticklabels(['10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{0}', '10^{1}', '10^{2}', '10^{3}'])
    # ax.tick_params(labelsize=14)
    # fig5.tight_layout()
    # fig5.savefig(config.decluster_data_fig_nna, dpi=300)
    #
    # # ----- Figure 4: 累积地震次数随时间变化 -----
    # # 原始目录
    # t_orig_list = []
    # len_data = data.shape[0]
    # for i in range(len_data):
    #     t_orig_list.append(datetime(int(data[i, 0]), int(data[i, 1]), int(data[i, 2]),
    #                                 int(data[i, 3]), int(data[i, 4]), int(data[i, 5])))
    # t_orig = mdates.date2num(t_orig_list)
    # num_eq = np.arange(1, len(t_orig) + 1)
    # num_eq_normalized = num_eq / num_eq[-1]
    # fig4, ax4 = plt.subplots(figsize=(8, 6))
    # ax4.plot(t_orig, num_eq_normalized, '-o', linewidth=1.5, markersize=2,
    #          color='black', label='Original Catalog')
    # # 去相关目录
    # t_decl_list = []
    # for i in range(Declustered_Catalog.shape[0]):
    #     t_decl_list.append(datetime(int(Declustered_Catalog[i, 0]), int(Declustered_Catalog[i, 1]),
    #                                 int(Declustered_Catalog[i, 2]), int(Declustered_Catalog[i, 3]),
    #                                 int(Declustered_Catalog[i, 4]), int(Declustered_Catalog[i, 5])))
    # t_decl = mdates.date2num(t_decl_list)
    # num_eq_decl = np.arange(1, len(t_decl) + 1)
    # num_eq_decl_normalized = num_eq_decl / num_eq_decl[-1] if len(num_eq_decl) > 0 else num_eq_decl
    # ax4.plot(t_decl, num_eq_decl_normalized, '-o', linewidth=1.5, markersize=2,
    #          color='red', label='Declustered Catalog')
    # ax4.set_xlabel('Time', fontsize=16, fontweight='bold')
    # ax4.set_ylabel('Normalized Cumulative Number of Earthquakes', fontsize=16, fontweight='bold')
    # ax4.set_title('Normalized Cumulative Earthquake Count Over Time', fontsize=18, fontweight='bold')
    # ax4.grid(True)
    # ax4.legend(fontsize=14, loc='upper left')
    # ax4.xaxis_date()
    # fig4.autofmt_xdate()
    # ax4.tick_params(labelsize=16)
    # fig4.tight_layout()
    # fig4.savefig(config.decluster_data_fig_cum, dpi=300)
    #
    # # ======== 将去相关目录输出到文件 ============
    # output_filename = config.data_select_decluster
    # with open(output_filename, 'w') as fid:
    #     # 格式为：'%d %02d %02d %02d %02d %05.2f %f %f %04.2f %04.2f %d\n'
    #     for i in range(Declustered_Catalog.shape[0]):
    #         yr = int(round(Declustered_Catalog[i, 0]))
    #         mon = int(round(Declustered_Catalog[i, 1]))
    #         day_ = int(round(Declustered_Catalog[i, 2]))
    #         hr = int(round(Declustered_Catalog[i, 3]))
    #         minu = int(round(Declustered_Catalog[i, 4]))
    #         sec_ = Declustered_Catalog[i, 5]
    #         lat_val = Declustered_Catalog[i, 6]
    #         lon_val = Declustered_Catalog[i, 7]
    #         dep_val = Declustered_Catalog[i, 8]
    #         mag_val = Declustered_Catalog[i, 9]
    #         extra = int(round(Declustered_Catalog[i, 10]))
    #         line = f"{yr} {mon:02d} {day_:02d} {hr:02d} {minu:02d} {sec_:05.2f} {lat_val:f} {lon_val:f} {dep_val:04.2f} {mag_val:04.2f} {extra:d}\n"
    #         fid.write(line)
