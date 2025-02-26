import numpy as np
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from scipy.ndimage import uniform_filter
from pyproj import Geod

#
# def haversine_distance(lat1, lon1, lat2, lon2, radius=6378.1):
#     """
#     计算大圆距离（单位：公里），支持标量或数组输入，角度单位为度。
#     """
#     # 转换为弧度
#     lat1 = np.radians(lat1)
#     lon1 = np.radians(lon1)
#     lat2 = np.radians(lat2)
#     lon2 = np.radians(lon2)
#     # 球面余弦公式
#     cos_val = np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2)
#     cos_val = np.clip(cos_val, -1, 1)
#     d = np.arccos(cos_val)
#     return d * radius

def distance(lat1, lon1, lat2, lon2):
    """
    Compute the great-circle distance (in km) and azimuths between two points
    given their latitude and longitude using the WGS84 ellipsoid.

    Parameters:
    lat1, lon1: Latitude and Longitude of the first point in degrees.
    lat2, lon2: Latitude and Longitude of the second point in degrees.

    Returns:
    dist_km: Distance between the two points in kilometers.
    azimuth1: Forward azimuth from point 1 to point 2 in degrees.
    azimuth2: Back azimuth from point 2 to point 1 in degrees.
    """
    geod = Geod(ellps="WGS84")
    azimuth1, azimuth2, dist_m = geod.inv(lon1, lat1, lon2, lat2)
    dist_km = dist_m / 1000.0  # Convert meters to kilometers
    return dist_km, azimuth1, azimuth2


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
        lat_scalar = np.full_like(lat[:ii], lat[ii])
        lon_scalar = np.full_like(lon[:ii], lon[ii])
        dr, az1, az2 = distance(lat_scalar, lon_scalar, lat[:ii], lon[:ii])
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