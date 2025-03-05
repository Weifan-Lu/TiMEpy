import matplotlib.pyplot as plt
import numpy as np
import main_function as mf
import visualization as vis
from matplotlib import rcParams

# rcParams['font.family'] = 'Arial'

def ana2_entire_region(config):
    # Close all existing figure windows
    plt.close('all')
    print('====== Analysis | Entire region: Start ======')
    print('Loading observational data (txt files)')

    # Read observational data (txt files)
    AllphaN = mf.load_txt_data(config.phase_stress_obs + '_N.txt')
    AllphaS = mf.load_txt_data(config.phase_stress_obs + '_S.txt')
    AllphaCFS = mf.load_txt_data(config.phase_stress_obs + '_CFS.txt')
    AllphaVol = mf.load_txt_data(config.phase_stress_obs + '_Vol.txt')

    print('Loading reference data (txt files)')

    # Read reference data (txt files)
    AllphaN_ref = mf.load_txt_data(config.phase_stress_ref + '_N.txt')
    AllphaS_ref = mf.load_txt_data(config.phase_stress_ref + '_S.txt')
    AllphaCFS_ref = mf.load_txt_data(config.phase_stress_ref + '_CFS.txt')
    AllphaVol_ref = mf.load_txt_data(config.phase_stress_ref + '_Vol.txt')

    # Use phase binning parameters from the configuration
    bin_phase = config.bin_phase

    # --- Step2. Instantaneous Phase Analysis ---
    # Assume the second column contains phase data (MATLAB AllphaX(:,2) → Python AllphaX[:,1])
    RPM_raVol, ph_shiftVol, _, _ = mf.modulation_phase(AllphaVol[:, 1] * 180 / np.pi, bin_phase,
                                                       AllphaVol_ref[:, 1] * 180 / np.pi, 'no')
    RPM_raS, ph_shiftS, _, _ = mf.modulation_phase(AllphaS[:, 1] * 180 / np.pi, bin_phase,
                                                     AllphaS_ref[:, 1] * 180 / np.pi, 'no')
    RPM_raN, ph_shiftN, _, _ = mf.modulation_phase(AllphaN[:, 1] * 180 / np.pi, bin_phase,
                                                     AllphaN_ref[:, 1] * 180 / np.pi, 'no')
    RPM_raCFS, ph_shiftCFS, _, _ = mf.modulation_phase(AllphaCFS[:, 1] * 180 / np.pi, bin_phase,
                                                       AllphaCFS_ref[:, 1] * 180 / np.pi, 'no')

    # --- Step3. Instantaneous Stress Modulation Analysis ---
    # Assume the third column contains stress data (MATLAB AllphaX(:,3) → Python AllphaX[:,2])
    initial_guess = config.initial_guess

    a_estimated_CFS, c_estimated_CFS, delta_a_CFS, delta_c_CFS, Stress_AM_bk_CFS, Stress_AM_CFS, shear_stress_CFS, shear_stress_kPa_CFS, event_rate_CFS = mf.modulation_stress(
        AllphaCFS_ref[:, 2], AllphaCFS[:, 2], initial_guess, 'no')

    a_estimated_S, c_estimated_S, delta_a_S, delta_c_S, Stress_AM_bk_S, Stress_AM_S, shear_stress_S, shear_stress_kPa_S, event_rate_S = mf.modulation_stress(
        AllphaS_ref[:, 2], AllphaS[:, 2], initial_guess, 'no')

    a_estimated_N, c_estimated_N, delta_a_N, delta_c_N, Stress_AM_bk_N, Stress_AM_N, shear_stress_N, shear_stress_kPa_N, event_rate_N = mf.modulation_stress(
        AllphaN_ref[:, 2], AllphaN[:, 2], initial_guess, 'no')

    a_estimated_Vol, c_estimated_Vol, delta_a_Vol, delta_c_Vol, Stress_AM_bk_Vol, Stress_AM_Vol, shear_stress_Vol, shear_stress_kPa_Vol, event_rate_Vol = mf.modulation_stress(
        AllphaVol_ref[:, 2], AllphaVol[:, 2], initial_guess, 'no')

    # --- Plotting and Saving ---
    # Figure 1: Volumetric stress phase modulation plot
    RPM_raVol, ph_shiftVol, _, _ = mf.modulation_phase(AllphaVol[:, 1] * 180 / np.pi, bin_phase,
                                                       AllphaVol_ref[:, 1] * 180 / np.pi, 'yes1')
    plt.savefig(config.TM_all_region_Vol, format='pdf')
    print(f"The tidal phase and amplitude (Vol) figure plotted and saved to {config.TM_all_region_Vol}")

    # Figure 2: Normal stress phase modulation plot
    RPM_raN, ph_shiftN, _, _ = mf.modulation_phase(AllphaN[:, 1] * 180 / np.pi, bin_phase,
                                                     AllphaN_ref[:, 1] * 180 / np.pi, 'yes1')
    plt.savefig(config.TM_all_region_N, format='pdf')
    print(f"The tidal phase and amplitude (Normal stress) figure plotted and saved to {config.TM_all_region_N}")

    # Figure 3: Shear stress phase modulation plot
    RPM_raS, ph_shiftS, _, _ = mf.modulation_phase(AllphaS[:, 1] * 180 / np.pi, bin_phase,
                                                     AllphaS_ref[:, 1] * 180 / np.pi, 'yes1')
    plt.savefig(config.TM_all_region_S, format='pdf')
    print(f"The tidal phase and amplitude (Shear stress) figure plotted and saved to {config.TM_all_region_S}")

    # Figure 4: Coulomb stress phase modulation plot
    RPM_raCFS, ph_shiftCFS, _, _ = mf.modulation_phase(AllphaCFS[:, 1] * 180 / np.pi, bin_phase,
                                                       AllphaCFS_ref[:, 1] * 180 / np.pi, 'yes1')
    plt.savefig(config.TM_all_region_CFS, format='pdf')
    print(f"The tidal phase and amplitude (Coulomb stress) figure plotted and saved to {config.TM_all_region_S}")

    fig = vis.plot_tidal_sensitivity_2x2(
        # Volumetric (First subplot)
        a_estimated_Vol, c_estimated_Vol, delta_a_Vol, delta_c_Vol,
        Stress_AM_bk_Vol, Stress_AM_Vol, shear_stress_Vol, shear_stress_kPa_Vol, event_rate_Vol,
        # Normal (Second subplot)
        a_estimated_N, c_estimated_N, delta_a_N, delta_c_N,
        Stress_AM_bk_N, Stress_AM_N, shear_stress_N, shear_stress_kPa_N, event_rate_N,
        # Shear (Third subplot)
        a_estimated_S, c_estimated_S, delta_a_S, delta_c_S,
        Stress_AM_bk_S, Stress_AM_S, shear_stress_S, shear_stress_kPa_S, event_rate_S,
        # Coulomb (Fourth subplot)
        a_estimated_CFS, c_estimated_CFS, delta_a_CFS, delta_c_CFS,
        Stress_AM_bk_CFS, Stress_AM_CFS, shear_stress_CFS, shear_stress_kPa_CFS, event_rate_CFS
    )
    # Save the figure to the path specified in the configuration file
    fig.savefig(config.TM_all_region_tidal_sensitivity_cfs, format='pdf')
    plt.close(fig)
    print(f"The tidal sensitivity figure plotted and saved to {config.TM_all_region_tidal_sensitivity_cfs}")




    # # Schuster 检验（这里只提供一个示例实现）
    # p_values_N = mf.schuster_test(AllphaN[:, 1])
    # p_values_S = mf.schuster_test(AllphaS[:, 1])
    # p_values_CFS = mf.schuster_test(AllphaCFS[:, 1])
    # p_values_Vol = mf.schuster_test(AllphaVol[:, 1])
    #
    # # Volumetric stress sensitivity 单独成图
    # plt.figure(figsize=(12, 12))
    # a_estimated_Vol, c_estimated_Vol, delta_a_Vol, delta_c_Vol= mf.modulation_stress(
    #     AllphaVol_ref[:, 2], AllphaVol[:, 2], initial_guess, 'yes')
    # plt.xlabel('Volumetric strain 10^{10}')
    # plt.savefig(config.TM_all_region_tidal_sensitivity_vol, format='pdf')
    # plt.close()
    #
    # # Normal stress sensitivity 单独成图
    # plt.figure(figsize=(12, 12))
    # a_estimated_N, c_estimated_N, delta_a_N, delta_c_N = mf.modulation_stress(
    #     AllphaN_ref[:, 2], AllphaN[:, 2], initial_guess, 'yes')
    # plt.xlabel('Normal stress (Pa)')
    # plt.savefig(config.TM_all_region_tidal_sensitivity_n, format='pdf')
    # plt.close()
    #
    # # Shear stress sensitivity 单独成图
    # plt.figure(figsize=(12, 12))
    # a_estimated_S, c_estimated_S, delta_a_S, delta_c_S = mf.modulation_stress(
    #     AllphaS_ref[:, 2], AllphaS[:, 2], initial_guess, 'yes')
    # plt.xlabel('Shear stress (Pa)')
    # plt.savefig(config.TM_all_region_tidal_sensitivity_s, format='pdf')
    # plt.close()
    #
    # # Coulomb stress sensitivity 单独成图
    # plt.figure(figsize=(12, 12))
    # a_estimated_CFS, c_estimated_CFS, delta_a_CFS, delta_c_CFS = mf.modulation_stress(
    #     AllphaCFS_ref[:, 2], AllphaCFS[:, 2], initial_guess, 'yes')
    # plt.xlabel('Coulomb stress (Pa)')
    # plt.savefig(config.TM_all_region_tidal_sensitivity_cfs, format='pdf')
    # plt.close()
