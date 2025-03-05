import numpy as np
import datetime
import matplotlib.dates as mdates
import main_function as mf
import visualization as vis

def ana1_calc_tidal_phase(config, opt):
    """
    Calculate tidal phases based on the parameter dictionary `params` and the mode `opt`, and save the results as a TXT file.

    The `params` dictionary must include the following entries (some parameters are provided in MATLAB datenum format, in units of days):
      - data_select_decluster : Earthquake catalog file (used if opt == '1')
      - data_select           : Earthquake catalog file (used if opt != '1')
      - stress_to_mat         : Path to the .mat file containing stress data (which includes the variables normal_stress, shear_stress, volumetric_strain)
      - t_stress_start        : Start time for the stress data (either a string in 'YYYY-MM-DD' format or a datenum value)
      - t_sample              : Sampling interval (in seconds)
      - miu                   : Friction coefficient
      - start_time            : Start time for the analysis (datenum value)
      - end_time              : End time for the analysis (datenum value)
      - t_search              : Search window (in days; note that the MATLAB code calculates the number of samples as t_search*86400/t_sample)
      - ref_interval          : Uniform reference time interval (in days)
      - phase_stress_obs      : Filename prefix for saving observation results as TXT files (the script automatically appends _N.txt, _S.txt, etc.)
      - phase_stress_ref      : Filename prefix for saving reference results as TXT files

    If opt is '1', data_select_decluster is used; otherwise, data_select is used.
    """
    ###
    print('====== Analysis | Earthquakes are associated with tidal stresses: Start ======')

    print('Start calculating the phase of the tide (Obs)')
    ###
    # --------------------
    # 1. Read earthquake catalog data
    if opt == '1':
        data = np.loadtxt(config.data_select_decluster)
    else:
        data = np.loadtxt(config.data_select)
    # Assume that each row in data has the following columns:
    # year, month, day, hour, minute, sec, (others), lat, lon, depth
    years = data[:, 0].astype(int)
    months = data[:, 1].astype(int)
    days = data[:, 2].astype(int)
    hours = data[:, 3].astype(int)
    minutes = data[:, 4].astype(int)
    secs = data[:, 5]
    # Construct a list of datetime objects
    eq_datetimes = [datetime.datetime(y, m, d, h, mi, int(s))
                    for y, m, d, h, mi, s in zip(years, months, days, hours, minutes, secs)]
    # Convert to datenum (matplotlib's date2num is largely compatible with MATLAB's datenum format)
    t_decluster = mdates.date2num(eq_datetimes)

    # --------------------
    # 2. Read stress data
    # Only read the required columns (assume the file has no header; if there is a header, add skiprows=1)
    data_stress = np.loadtxt(config.output_stress_txt, usecols=(3, 4, 5))

    # According to the file format, the indices correspond to:
    volumetric_strain = data_stress[:, 0]
    shear_stress = data_stress[:, 1]
    normal_stress = data_stress[:, 2]
    miu = config.miu
    StressN = normal_stress
    StressS = shear_stress
    StressVol = volumetric_strain
    StressCFS = shear_stress + miu * normal_stress

    # --------------------
    # 3. Construct the stress data time series ttide
    # t_stress_start can be a string or a numeric value
    if isinstance(config.t_stress_start, str):
        t0_dt = datetime.datetime.strptime(config.t_stress_start, '%Y-%m-%d')
        t_stress_start = mdates.date2num(t0_dt)
    else:
        t_stress_start = config.t_stress_start
    t_sample = config.t_sample  # Unit: seconds
    n_stress = len(normal_stress)
    ttide = t_stress_start + np.arange(n_stress) * t_sample / 86400.0  # Convert to days
    t_tide = ttide.copy()

    # --------------------
    # 4. Calculate tidal phases for catalog data (obs)
    t_start = config.start_time  # datenum value
    t_end = config.end_time
    print('start time:', t_start, '  end time: ', t_end, ' len_stress: ', len(StressN), 'len_time: ', len(ttide))
    # If t_start and t_end are datetime objects, convert them to datenum format (float)
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

    for ix in indices:
        t_cut = t_decluster[ix]
        Index = np.argmin(np.abs(t_tide - t_cut))
        range_indices = np.arange(Index - incat, Index + incat + 1)
        # print(range_indices)
        # Check if the index range is out of bounds
        if range_indices[0] < 0 or range_indices[-1] >= len(StressN):
            # print('wrong')
            continue
        StressN_cut = StressN[range_indices]
        StressS_cut = StressS[range_indices]
        StressVol_cut = StressVol[range_indices]
        StressCFS_cut = StressCFS[range_indices]
        ttide_cut = t_tide[range_indices]

        # Calculate tidal phase and other information separately (here only using phase, stress value, and stress rate)
        phaN, _, _, _, stress_at_t_N, stress_rate_at_t_N = mf.calculate_tidal_phase(ttide_cut, StressN_cut, t_cut)
        phaS, _, _, _, stress_at_t_S, stress_rate_at_t_S = mf.calculate_tidal_phase(ttide_cut, StressS_cut, t_cut)
        phaCFS, _, _, _, stress_at_t_CFS, stress_rate_at_t_CFS = mf.calculate_tidal_phase(ttide_cut, StressCFS_cut,
                                                                                          t_cut)
        phaVol, _, _, _, stress_at_t_Vol, stress_rate_at_t_Vol = mf.calculate_tidal_phase(ttide_cut, StressVol_cut,
                                                                                          t_cut)
        dt = mdates.num2date(t_cut)
        lon_val = data[ix, 7]
        lat_val = data[ix, 6]
        depth_val = data[ix, 8]
        mag_val = data[ix, 9]
        obs_list_N.append(
            [t_cut, phaN, stress_at_t_N, stress_rate_at_t_N, dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second,
             lat_val, lon_val, depth_val, mag_val])
        obs_list_S.append(
            [t_cut, phaS, stress_at_t_S, stress_rate_at_t_S, dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second,
             lat_val, lon_val, depth_val, mag_val])
        obs_list_CFS.append(
            [t_cut, phaCFS, stress_at_t_CFS, stress_rate_at_t_CFS, dt.year, dt.month, dt.day, dt.hour, dt.minute,
             dt.second, lat_val, lon_val, depth_val, mag_val])
        obs_list_Vol.append(
            [t_cut, phaVol, stress_at_t_Vol, stress_rate_at_t_Vol, dt.year, dt.month, dt.day, dt.hour, dt.minute,
             dt.second, lat_val, lon_val, depth_val, mag_val])
    AllphaN = np.array(obs_list_N)
    AllphaS = np.array(obs_list_S)
    AllphaCFS = np.array(obs_list_CFS)
    AllphaVol = np.array(obs_list_Vol)

    print('Plot stress v.s. eqs')

    vis.plot_stress_earthquake(AllphaVol[:, 0], AllphaVol[:, 2], AllphaN[:, 2], AllphaS[:, 2], AllphaCFS[:, 2], t_tide,
                               StressVol, StressN, StressS, StressCFS, config)

    print(f"The figure plotted and saved to {config.TM_all_region_stress}")

    print('Start calculating the phase of the tide (Ref)')

    # --------------------
    # 5. Generate a uniform reference time series and calculate tidal phases (ref)
    ref_interval = config.ref_interval  # Unit: days
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
        phaCFS, _, _, _, stress_at_t_CFS, stress_rate_at_t_CFS = mf.calculate_tidal_phase(ttide_cut, StressCFS_cut,
                                                                                          t_cut)
        phaVol, _, _, _, stress_at_t_Vol, stress_rate_at_t_Vol = mf.calculate_tidal_phase(ttide_cut, StressVol_cut,
                                                                                          t_cut)

        ref_list_N.append([t_cut, phaN, stress_at_t_N, stress_rate_at_t_N])
        ref_list_S.append([t_cut, phaS, stress_at_t_S, stress_rate_at_t_S])
        ref_list_CFS.append([t_cut, phaCFS, stress_at_t_CFS, stress_rate_at_t_CFS])
        ref_list_Vol.append([t_cut, phaVol, stress_at_t_Vol, stress_rate_at_t_Vol])
    AllphaN_ref = np.array(ref_list_N)
    AllphaS_ref = np.array(ref_list_S)
    AllphaCFS_ref = np.array(ref_list_CFS)
    AllphaVol_ref = np.array(ref_list_Vol)

    # --------------------
    # 6. Save the results as TXT files
    header = "t (datenum)\tphase (rad)\tstress_at_t\tstress_rate_at_t\tyear\tmonth\tday\thour\tmin\tsec\tlat\tlon\tdep\tmag"
    np.savetxt(config.phase_stress_obs + '_N.txt', AllphaN, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_obs + '_S.txt', AllphaS, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_obs + '_CFS.txt', AllphaCFS, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_obs + '_Vol.txt', AllphaVol, fmt='%.8f', delimiter='\t', header=header)
    #
    np.savetxt(config.phase_stress_ref + '_N.txt', AllphaN_ref, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_ref + '_S.txt', AllphaS_ref, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_ref + '_CFS.txt', AllphaCFS_ref, fmt='%.8f', delimiter='\t', header=header)
    np.savetxt(config.phase_stress_ref + '_Vol.txt', AllphaVol_ref, fmt='%.8f', delimiter='\t', header=header)

    print("Tidal phase results have been saved as TXT files with the following details:")
    print(f"  - Observation files:")
    print(f"      {config.phase_stress_obs + '_N.txt'} (Normal Stress)")
    print(f"      {config.phase_stress_obs + '_S.txt'} (Shear Stress)")
    print(f"      {config.phase_stress_obs + '_CFS.txt'} (Coulomb Failure Stress)")
    print(f"      {config.phase_stress_obs + '_Vol.txt'} (Volume Change)")
    print(f"  - Reference files:")
    print(f"      {config.phase_stress_ref + '_N.txt'} (Normal Stress)")
    print(f"      {config.phase_stress_ref + '_S.txt'} (Shear Stress)")
    print(f"      {config.phase_stress_ref + '_CFS.txt'} (Coulomb Failure Stress)")
    print(f"      {config.phase_stress_ref + '_Vol.txt'} (Volume Change)")

    print('====== Analysis | Earthquakes are associated with tidal stresses: End ======')
