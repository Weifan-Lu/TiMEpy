import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import main_function as mf
import visualization as vis
import os


def ana3_temp_variation(config):
    plt.close('all')

    # Read observational data (txt files)
    AllphaN = mf.load_txt_data(config.phase_stress_obs + '_N.txt')
    AllphaS = mf.load_txt_data(config.phase_stress_obs + '_S.txt')
    AllphaCFS = mf.load_txt_data(config.phase_stress_obs + '_CFS.txt')
    AllphaVol = mf.load_txt_data(config.phase_stress_obs + '_Vol.txt')

    # Read reference data (txt files)
    AllphaN_ref = mf.load_txt_data(config.phase_stress_ref + '_N.txt')
    AllphaS_ref = mf.load_txt_data(config.phase_stress_ref + '_S.txt')
    AllphaCFS_ref = mf.load_txt_data(config.phase_stress_ref + '_CFS.txt')
    AllphaVol_ref = mf.load_txt_data(config.phase_stress_ref + '_Vol.txt')

    # ----------------------------
    # 2. Get time range and parameter settings
    dt_start = config.start_time
    dt_end = config.end_time
    tmin = mdates.date2num(dt_start)
    tmax = mdates.date2num(dt_end)
    step_days = config.step_days
    window_days = config.window_days
    print(tmin,tmax)

    # Initialize various result lists
    p_values_N = []
    p_values_S = []
    p_values_CFS = []
    p_values_Vol = []
    Amp_S = []
    Amp_N = []
    Amp_Vol = []
    Amp_CFS = []
    shift_S = []
    shift_N = []
    shift_Vol = []
    shift_CFS = []
    A_est_S = []
    A_est_CFS = []
    A_est_N = []
    A_est_Vol = []
    time_series = []

    # Here, `all_syn_hist` in MATLAB is pre-allocated as a cell array. In Python, a list can be used directly.
    all_syn_hist = []
    bin_phase = config.bin_phase
    atimes = 3000
    idxc = 0

    # ----------------------------
    # 3. Sliding window processing
    # Loop from tmin to tmax - window_days (assuming `datenum` format is a floating-point number)
    for t_start in np.arange(tmin, tmax - window_days + 1e-6, step_days):
        t_end = t_start + window_days

        # Select observational data within the current window
        indices = np.where((AllphaS[:, 0] >= t_start) & (AllphaS[:, 0] <= t_end))[0]
        indices_ref = np.where((AllphaS_ref[:, 0] >= t_start) & (AllphaS_ref[:, 0] <= t_end))[0]

        # Calculate the center time of the window (here, we directly take the window endpoint; it can be modified to t_start + window_days/2 if needed)
        t_cent = t_start + window_days
        time_series.append(t_cent)
        idxc += 1
        number_eq = len(indices)

        # Compute modulation parameters (Note: MATLAB indexing starts from 1, whereas Python indexing starts from 0)
        # Assuming the second column (index 1) of `Allpha*` is the phase (unit: radians), convert it to degrees before passing in.

        RPM_raS, ph_shiftS, _, _  = mf.modulation_phase(AllphaS[indices, 1] * 180 / np.pi,
                                             bin_phase,
                                             AllphaS_ref[indices_ref, 1] * 180 / np.pi,
                                             'no')
        RPM_raN, ph_shiftN, _, _ = mf.modulation_phase(AllphaN[indices, 1] * 180 / np.pi,
                                             bin_phase,
                                             AllphaN_ref[indices_ref, 1] * 180 / np.pi,
                                             'no')
        RPM_raCFS, ph_shiftCFS, _, _ = mf.modulation_phase(AllphaCFS[indices, 1] * 180 / np.pi,
                                                 bin_phase,
                                                 AllphaCFS_ref[indices_ref, 1] * 180 / np.pi,
                                                 'no')
        RPM_raVol, ph_shiftVol, _, _ = mf.modulation_phase(AllphaVol[indices, 1] * 180 / np.pi,
                                                 bin_phase,
                                                 AllphaVol_ref[indices_ref, 1] * 180 / np.pi,
                                                 'no')

        # If some indices are large, generate detailed figures (mode "yes3")
        if RPM_raCFS > 1.11:
            output_dir = config.output_fig
            # Figure (1): Volumetric strain
            plt.figure(1, figsize=(60, 16))
            RPM_raVol, ph_shiftVol, _, _  = mf.modulation_phase(AllphaVol[indices, 1] * 180 / np.pi,bin_phase,AllphaVol_ref[indices_ref, 1] * 180 / np.pi,'yes3')
            # Convert `t_start` (datenum) to a datetime object
            date_time = mdates.num2date(t_start)
            filename = os.path.join(output_dir, f"Result_tidal_amp_phase_Vol_{date_time.strftime('%Y%m%d%H%M%S')}.pdf")
            plt.savefig(filename, format='pdf')

            # figure (2): Normal stress
            plt.figure(2, figsize=(60, 16))
            RPM_raN, ph_shiftN , _, _= mf.modulation_phase(AllphaN[indices, 1] * 180 / np.pi,bin_phase,AllphaN_ref[indices_ref, 1] * 180 / np.pi,'yes3')
            date_time = mdates.num2date(t_start)
            filename = os.path.join(output_dir, f"Result_tidal_amp_phase_N_{date_time.strftime('%Y%m%d%H%M%S')}.pdf")
            plt.savefig(filename, format='pdf')

            # figure (3): Shear stress
            plt.figure(3, figsize=(60, 16))
            RPM_raS, ph_shiftS , _, _= mf.modulation_phase(AllphaS[indices, 1] * 180 / np.pi,bin_phase,AllphaS_ref[indices_ref, 1] * 180 / np.pi,'yes3')
            date_time = mdates.num2date(t_start)
            filename = os.path.join(output_dir, f"Result_tidal_amp_phase_S_{date_time.strftime('%Y%m%d%H%M%S')}.pdf")
            plt.savefig(filename, format='pdf')

            # figure (4): Coulomb stress
            plt.figure(4, figsize=(60, 16))
            RPM_raCFS, ph_shiftCFS, _, _= mf.modulation_phase(AllphaCFS[indices, 1] * 180 / np.pi,bin_phase,AllphaCFS_ref[indices_ref, 1] * 180 / np.pi,'yes1')
            date_time = mdates.num2date(t_start)
            filename = os.path.join(output_dir, f"Result_tidal_amp_phase_CFS_{date_time.strftime('%Y%m%d%H%M%S')}.pdf")
            plt.savefig(filename, format='pdf')
            plt.close('all')

            # Coulomb stress sensitivity as a separate figure
            plt.figure(figsize=(6, 6))
            a_estimated_CFS, c_estimated_CFS, delta_a_CFS, delta_c_CFS, Stress_AM_bk_CFS, Stress_AM_CFS, shear_stress_CFS, shear_stress_kPa_CFS, event_rate_CFS = mf.modulation_stress(
                AllphaCFS_ref[indices_ref, 2], AllphaCFS[indices, 2], initial_guess, 'yes')
            plt.xlabel('Coulomb stress (Pa)')
            filename = os.path.join(output_dir, f"Result_tidal_sensitivity_CFS_{date_time.strftime('%Y%m%d%H%M%S')}.pdf")
            plt.savefig(filename, format='pdf')

            plt.close()

        # Compute p-values (using Schuster test function, input is phase data in radians)
        p_values_N.append(mf.schuster_test(AllphaN[indices, 1]))
        p_values_S.append(mf.schuster_test(AllphaS[indices, 1]))
        p_values_CFS.append(mf.schuster_test(AllphaCFS[indices, 1]))
        p_values_Vol.append(mf.schuster_test(AllphaVol[indices, 1]))

        # Save amplitude and phase shift information
        Amp_S.append(RPM_raS)
        Amp_N.append(RPM_raN)
        Amp_Vol.append(RPM_raVol)
        Amp_CFS.append(RPM_raCFS)
        shift_S.append(ph_shiftS)
        shift_N.append(ph_shiftN)
        shift_Vol.append(ph_shiftVol)
        shift_CFS.append(ph_shiftCFS)

        # Modulation stress fitting, the 3rd column in MATLAB corresponds to stress values,
        # which corresponds to index 2 in Python (0-based)
        initial_guess = [2, 0.1]
        a_estimated_S, c_estimated_S, delta_a_S, delta_c_S, Stress_AM_bk_S, Stress_AM_S, shear_stress_S, shear_stress_kPa_S, event_rate_S = mf.modulation_stress(
            AllphaS_ref[indices_ref, 2], AllphaS[indices, 2], initial_guess, 'no')
        a_estimated_N, c_estimated_N, delta_a_N, delta_c_N, Stress_AM_bk_N, Stress_AM_N, shear_stress_N, shear_stress_kPa_N, event_rate_N = mf.modulation_stress(
            AllphaN_ref[indices_ref, 2], AllphaN[indices, 2], initial_guess, 'no')
        a_estimated_Vol, c_estimated_Vol, delta_a_Vol, delta_c_Vol, Stress_AM_bk_Vol, Stress_AM_Vol, shear_stress_Vol, shear_stress_kPa_Vol, event_rate_Vol = mf.modulation_stress(
            AllphaVol_ref[indices_ref, 2], AllphaVol[indices, 2], initial_guess, 'no')
        a_estimated_CFS, c_estimated_CFS, delta_a_CFS, delta_c_CFS, Stress_AM_bk_CFS, Stress_AM_CFS, shear_stress_CFS, shear_stress_kPa_CFS, event_rate_CFS = mf.modulation_stress(
            AllphaCFS_ref[indices_ref, 2], AllphaCFS[indices, 2], initial_guess, 'no')

        A_est_S.append([a_estimated_S, c_estimated_S, delta_a_S, delta_c_S, number_eq])
        A_est_N.append([a_estimated_N, c_estimated_N, delta_a_N, delta_c_N, number_eq])
        A_est_Vol.append([a_estimated_Vol, c_estimated_Vol, delta_a_Vol, delta_c_Vol, number_eq])
        A_est_CFS.append([a_estimated_CFS, c_estimated_CFS, delta_a_CFS, delta_c_CFS, number_eq])

    # In MATLAB code, empty rows are removed; here, all_syn_hist can be processed as needed

    # ----------------------------
    # 4. Plot results
    komo_t = config.mainshock_t

    # Precompute x-axis year labels (time_series stores datenum values)
    xlim_range = [min(time_series), max(time_series)]
    start_year = mdates.num2date(xlim_range[0]).year
    end_year = mdates.num2date(xlim_range[1]).year
    years = np.arange(start_year, end_year + 1)

    # (1) p-value plots — four subplots
    fig, axs = plt.subplots(4, 1, figsize=(8, 12))
    vis.plot_p_value_subplot(0, p_values_Vol, time_series, komo_t, 'Volumetric strain', ax=axs[0])
    vis.plot_p_value_subplot(1, p_values_N, time_series, komo_t, 'Normal stress', ax=axs[1])
    vis.plot_p_value_subplot(2, p_values_S, time_series, komo_t, 'Shear stress', ax=axs[2])
    vis.plot_p_value_subplot(3, p_values_CFS, time_series, komo_t, 'Coulomb stress', ax=axs[3])
    plt.tight_layout()
    plt.savefig(config.TM_time_tidal_P_value, format='pdf')

    # (2) Amplitude plots — four subplots
    fig, axs = plt.subplots(4, 1, figsize=(8, 12))
    vis.plot_amp_value_subplot(0, Amp_Vol, time_series, komo_t, 'Volumetric strain', ax=axs[0])
    vis.plot_amp_value_subplot(1, Amp_N, time_series, komo_t, 'Normal stress', ax=axs[1])
    vis.plot_amp_value_subplot(2, Amp_S, time_series, komo_t, 'Shear stress', ax=axs[2])
    vis.plot_amp_value_subplot(3, Amp_CFS, time_series, komo_t, 'Coulomb stress', ax=axs[3])
    plt.tight_layout()
    plt.savefig(config.TM_time_tidal_amplitude, format='pdf')

    # (3) Phase shift plots — four subplots
    fig, axs = plt.subplots(4, 1, figsize=(8, 12))
    vis.plot_pha_value_subplot(0, shift_Vol, time_series, komo_t, 'Volumetric strain', ax=axs[0])
    vis.plot_pha_value_subplot(1, shift_N, time_series, komo_t, 'Normal stress', ax=axs[1])
    vis.plot_pha_value_subplot(2, shift_S, time_series, komo_t, 'Shear stress', ax=axs[2])
    vis.plot_pha_value_subplot(3, shift_CFS, time_series, komo_t, 'Coulomb stress', ax=axs[3])
    plt.tight_layout()
    plt.savefig(config.TM_time_tidal_phase_shift, format='pdf')

    # (4) Sensitivity plots — four subplots
    fig, axs = plt.subplots(4, 1, figsize=(8, 12))
    vis.plot_a_value_subplot(0, A_est_Vol, time_series, komo_t, 'Volumetric strain', ax=axs[0])
    vis.plot_a_value_subplot(1, A_est_N, time_series, komo_t, 'Normal stress', ax=axs[1])
    vis.plot_a_value_subplot(2, A_est_S, time_series, komo_t, 'Shear stress', ax=axs[2])
    vis.plot_a_value_subplot(3, A_est_CFS, time_series, komo_t, 'Coulomb stress', ax=axs[3])
    plt.tight_layout()
    plt.savefig(config.TM_time_tidal_sensitivity, format='pdf')

    # (5) Comprehensive plots — four subplots (example with Coulomb stress)
    fig, axs = plt.subplots(4, 1, figsize=(8, 12))
    vis.plot_p_value_subplot(0, p_values_CFS, time_series, komo_t, '', ax=axs[0])
    vis.plot_pha_value_subplot(1, shift_CFS, time_series, komo_t, '', ax=axs[1])
    vis.plot_amp_value_subplot(2, Amp_CFS, time_series, komo_t, '', ax=axs[2])
    # vis.plot_pha_value_subplot(2, shift_CFS, time_series, komo_t, 'Coulomb stress', ax=axs[2])
    vis.plot_a_value_subplot(3, A_est_CFS, time_series, komo_t, 'Coulomb stress', ax=axs[3])
    plt.tight_layout()
    plt.savefig(config.TM_time_tidal_all_CFS, format='pdf')

    # (6) Earthquake event count plot (taking the 5th column of A_est_N, which is the event count)
    fig, ax = plt.subplots(1, 1, figsize=(20, 3))
    # Extract event count (Note: A_est_N is a list where each entry is [a_est, c_est, delta_a, delta_c, number_eq])
    number_events = [entry[4] for entry in A_est_N]
    ax.plot(time_series, number_events, '-o', linewidth=1.5, markersize=6, markerfacecolor='b')
    ax.set_xlabel('Time', fontsize=14)
    ax.set_ylabel('Number of Events', fontsize=14)
    ax.grid(True)
    ax.tick_params(axis='both', labelsize=18)
    # Use matplotlib.dates to format x-axis as years
    ax.xaxis_date()
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    plt.tight_layout()
    plt.savefig(config.TM_time_earthquake_number, format='pdf')

    plt.close('all')