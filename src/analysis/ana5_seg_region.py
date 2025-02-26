import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np
import main_function as mf
import visualization as vis
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'


def ana5_seg_region(config):
    # 关闭所有已有图窗
    plt.close('all')

    # 读取观测数据（txt 文件）
    AllphaN = mf.load_txt_data(config.phase_stress_obs + '_N.txt')
    AllphaS = mf.load_txt_data(config.phase_stress_obs + '_S.txt')
    AllphaCFS = mf.load_txt_data(config.phase_stress_obs + '_CFS.txt')
    AllphaVol = mf.load_txt_data(config.phase_stress_obs + '_Vol.txt')

    # 读取参考数据（txt 文件）
    AllphaN_ref = mf.load_txt_data(config.phase_stress_ref + '_N.txt')
    AllphaS_ref = mf.load_txt_data(config.phase_stress_ref + '_S.txt')
    AllphaCFS_ref = mf.load_txt_data(config.phase_stress_ref + '_CFS.txt')
    AllphaVol_ref = mf.load_txt_data(config.phase_stress_ref + '_Vol.txt')
    obs_vars = [AllphaN, AllphaS, AllphaCFS, AllphaVol]
    ref_vars = [AllphaN_ref, AllphaS_ref, AllphaCFS_ref, AllphaVol_ref]
    data = np.loadtxt(config.data_select_decluster)
    # Define stress types and extract corresponding variables from the loaded data.
    stress_types = ['Normal Stress', 'Shear Stress', 'CFS', 'Volumetric Strain']
    # Parameter settings
    grid_spacing = 0.2  # 每个网格 0.2 度
    threshold = 20  # 每个网格至少需要 20 个样本 (MATLAB comment said 10, but variable is 20)
    bin_phase = config.bin_phase
    initial_guess = config.initial_guess

    # Loop over each stress type and generate four subplots per figure.
    for i, stress_type in enumerate(stress_types):
        # Compute grid centers for the current stress type using a custom function.
        grid_centers = mf.grid_modulation_results_2D(
            obs_vars[i],
            ref_vars[i],
            bin_phase,
            initial_guess,
            grid_spacing,
            threshold
        )

        # Create a new figure with the desired size.
        plt.figure(figsize=(18, 12))

        # Plot first subplot: p_value
        ax = plt.subplot(2, 2, 1)
        vis.plot_filled_grids_p_value(grid_centers, grid_spacing, stress_type, data,ax)

        # Plot second subplot: a_value
        ax = plt.subplot(2, 2, 2)
        vis.plot_filled_grids_a_value(grid_centers, grid_spacing, stress_type, data,ax)

        # Plot third subplot: amp_value
        ax = plt.subplot(2, 2, 3)
        vis.plot_filled_grids_amp_value(grid_centers, grid_spacing, stress_type, data,ax)

        # Plot fourth subplot: pha_value
        ax = plt.subplot(2, 2, 4)
        vis.plot_filled_grids_pha_value(grid_centers, grid_spacing, stress_type, data,ax)

        # Optionally adjust the layout and display the figure.
        plt.tight_layout()
        plt.show()