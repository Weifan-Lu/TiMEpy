import numpy as np
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from scipy.ndimage import uniform_filter
from pyproj import Geod


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
    Perform a moving average smoothing on the 2D array A with a window size of (Nr, Nc).
    """
    return uniform_filter(A, size=(Nr, Nc))


def decluster_nna(start_time,data,df_val,b,p,q,eta0):
    """
        Process the earthquake catalog data based on control parameters. This function encapsulates
        the computation of Nij, Tij, Rij, the construction of BgLink/StLink, and family identification.

        Parameters:
            params (dict): A dictionary of control parameters that must include the following keys:
                - 'start_time': (optional) Start time (not used in this example).
                - 'data_select': Path to the data file (used when `data` is not provided).
                - 'df': Distance power exponent (e.g., 1.6).
                - 'b': b-value parameter.
                - 'p': p-value parameter.
                - 'q': q-value parameter.
                - 'eta0': Threshold for classification (compared with the logarithm of the third column of LinkDis).
            data (np.ndarray, optional): Data array. If provided, `data_select` is ignored.
                The data must have a shape of (N, 10) with columns in the following order:
                [year, month, day, hour, minute, second, lat, lon, depth, magnitude].

        Returns:
            dict: A dictionary containing the processed results with the following keys:
                - 'Declustered_Catalog': The declustered catalog. Each row corresponds to the original
                  `data` row, with an additional column for t/365 (in years).
                - 'Dataobj': The processed Dataobj array (first and second columns are logarithmic,
                  the third column is min_nj).
                - 'LinkDis': The LinkDis array (third column is logarithmic).
                - 'BgLink': The BgLink array.
                - 'StLink': The StLink array.
                - 'family_weak': The weak family array.
                - 'Mainshock': A dictionary of mainshock families, where each key corresponds to an
                  array of event indices in the family.
                - 'All_Mainshock': A list of all mainshock event indices.
                - 'All_Foreshock': A list of foreshock groups (each group is a 2D array).
                - 'All_Aftershock': A list of aftershock groups (each group is a 2D array).
                - 'All_Family': A merged array of family event indices.
                - 'All_Single': A set of single (isolated) events.
                - 'z': The third column of Dataobj (min_nj).
        """

    year = data[:, 0].astype(int)
    month = data[:, 1].astype(int)
    day = data[:, 2].astype(int)
    hour = data[:, 3].astype(int)
    minute = data[:, 4].astype(int)
    sec = data[:, 5]
    lat = data[:, 6]
    lon = data[:, 7]
    dep = data[:, 8]
    mag = data[:, 9]

    # Convert time information to a datetime object, then convert it to a matplotlib date number (unit: days)
    t_list = [datetime(y, m, d, h, mi, int(s))
              for y, m, d, h, mi, s in zip(year, month, day, hour, minute, sec)]
    t = mdates.date2num(t_list)

    # 数据行数
    len_data = data.shape[0]

    # ======== Calculate Nij, Tij, Rij, Dataobj, LinkDis ========
    Dataobj = []  # Store [min_tj, min_rj, min_nj]
    LinkDis = []  # Store [current event index, nearest neighbor event index, min_nj]
    for ii in range(1, len_data):  # ii from 1 to len_data-1
        dt = (t[ii] - t[:ii]) / 365.0  # Time difference (years)
        lat_scalar = np.full_like(lat[:ii], lat[ii])
        lon_scalar = np.full_like(lon[:ii], lon[ii])
        dr, az1, az2 = distance(lat_scalar, lon_scalar, lat[:ii], lon[:ii])
        dm = mag[:ii]
        Nij = dt * (dr ** df_val) * (10 ** (-b * dm))
        Tij = dt * (10 ** (-q * b * dm))
        Rij = (dr ** df_val) * (10 ** (-p * b * dm))
        # Replace elements with a value of 0 with nan
        Nij = np.where(Nij == 0, np.nan, Nij)
        Tij = np.where(Tij == 0, np.nan, Tij)
        Rij = np.where(Rij == 0, np.nan, Rij)
        if np.all(np.isnan(Nij)):
            continue
        # Find the minimum value and its index in Nij (ignoring nan)
        xnumx = np.nanargmin(Nij)
        min_nj = Nij[xnumx]
        min_tj = Tij[xnumx]
        min_rj = Rij[xnumx]
        Dataobj.append([min_tj, min_rj, min_nj])
        LinkDis.append([ii, xnumx, min_nj])
    Dataobj = np.array(Dataobj)
    LinkDis = np.array(LinkDis)

    # Take the logarithm of the third column in LinkDis
    if LinkDis.size > 0:
        LinkDis[:, 2] = np.log10(LinkDis[:, 2])

    # ======== Construct BgLink and StLink ============
    BgLink = []
    StLink = []
    for i in range(LinkDis.shape[0]):
        # If log10(min_nj) > eta0, it is classified as BgLink, otherwise it is classified as StLink
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

    # ======== family_weak ============
    family_weak = []
    if BgLink.size > 0 and StLink.size > 0:
        for row in BgLink:
            if int(row[0]) in StLink[:, 1].astype(int):
                family_weak.append([int(row[0]), 1])
    family_weak = np.array(family_weak) if len(family_weak) > 0 else np.empty((0, 2))

    # ======== Construct family index =============
    num_st = StLink.shape[0]
    Family_idx = np.arange(num_st)
    for i in range(1, num_st):
        for j in range(i):
            if (StLink[j, 0] == StLink[i, 1]) or (StLink[j, 1] == StLink[i, 1]):
                Family_idx[i] = Family_idx[j]
                break
    Fm = np.unique(Family_idx)

    # ======== Construct Mainshock  ============
    Mainshock = {}
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
        # Merge all event indexes in the family and remove duplicates
        AA = np.unique(np.concatenate((aa[:, 0], aa[:, 1])).astype(int))
        # Select the event with the largest magnitude as the main shock
        mags = mag[AA]
        num_max = np.argmax(mags)
        mainshock_idx = AA[num_max]
        All_Mainshock.append(mainshock_idx)
        # Construct a two-dimensional array, the second column is the main shock index
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

    # Merge event indexes in each family
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

    # ======== Single event =============
    # Note: 1 in MATLAB corresponds to 0 in Python (first event)
    single_candidates = set(BgLink[:, 0].astype(int)) if BgLink.size > 0 else set()
    weak_family = set(family_weak[:, 0].astype(int)) if family_weak.size > 0 else set()
    All_Single = {0} | (single_candidates - weak_family)

    # ======== Go to related events =============
    All_decluster = np.unique(np.concatenate((np.array(list(All_Single)), np.array(All_Mainshock))))
    Declustered_Catalog = []
    for idx in All_decluster:
        Declustered_Catalog.append(np.concatenate((data[int(idx), :], [t[int(idx)] / 365.0])))
    Declustered_Catalog = np.array(Declustered_Catalog)

    # ======== Dataobj  ============
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