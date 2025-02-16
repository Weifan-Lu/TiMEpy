import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np
import main_function as mf
import visualization as vis
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'


def ana2_all_region(config):
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

    # 使用配置中的相位分箱参数
    bin_phase = config.bin_phase

    # --- Step2. 瞬时相位分析 ---
    # 假设第二列为相位数据（MATLAB AllphaX(:,2) → Python AllphaX[:,1]）
    RPM_raVol, ph_shiftVol, _, _ = mf.modulation_phase(AllphaVol[:, 1] * 180 / np.pi, bin_phase,
                                                       AllphaVol_ref[:, 1] * 180 / np.pi, 'no')
    RPM_raS, ph_shiftS, _, _ = mf.modulation_phase(AllphaS[:, 1] * 180 / np.pi, bin_phase,
                                                     AllphaS_ref[:, 1] * 180 / np.pi, 'no')
    RPM_raN, ph_shiftN, _, _ = mf.modulation_phase(AllphaN[:, 1] * 180 / np.pi, bin_phase,
                                                     AllphaN_ref[:, 1] * 180 / np.pi, 'no')
    RPM_raCFS, ph_shiftCFS, _, _ = mf.modulation_phase(AllphaCFS[:, 1] * 180 / np.pi, bin_phase,
                                                       AllphaCFS_ref[:, 1] * 180 / np.pi, 'no')


    # Schuster 检验（这里只提供一个示例实现）
    p_values_N = mf.schuster_test(AllphaN[:, 1])
    p_values_S = mf.schuster_test(AllphaS[:, 1])
    p_values_CFS = mf.schuster_test(AllphaCFS[:, 1])
    p_values_Vol = mf.schuster_test(AllphaVol[:, 1])

    # --- Step3. 瞬时应力调制分析 ---
    # 假设第三列为应力数据（MATLAB AllphaX(:,3) → Python AllphaX[:,2]）
    initial_guess = config.initial_guess
    # 以下代码如有需要可启用：
    a_estimated_CFS, c_estimated_CFS, delta_a_CFS, delta_c_CFS, Stress_AM_bk_CFS, Stress_AM_CFS, shear_stress_CFS, shear_stress_kPa_CFS, event_rate_CFS = mf.modulation_stress(
        AllphaCFS_ref[:, 2], AllphaCFS[:, 2], initial_guess, 'no')

    a_estimated_S, c_estimated_S, delta_a_S, delta_c_S, Stress_AM_bk_S, Stress_AM_S, shear_stress_S, shear_stress_kPa_S, event_rate_S = mf.modulation_stress(
        AllphaS_ref[:, 2], AllphaS[:, 2], initial_guess, 'no')

    a_estimated_N, c_estimated_N, delta_a_N, delta_c_N, Stress_AM_bk_N, Stress_AM_N, shear_stress_N, shear_stress_kPa_N, event_rate_N = mf.modulation_stress(
        AllphaN_ref[:, 2], AllphaN[:, 2], initial_guess, 'no')

    a_estimated_Vol, c_estimated_Vol, delta_a_Vol, delta_c_Vol, Stress_AM_bk_Vol, Stress_AM_Vol, shear_stress_Vol, shear_stress_kPa_Vol, event_rate_Vol = mf.modulation_stress(
        AllphaVol_ref[:, 2], AllphaVol[:, 2], initial_guess, 'no')

    # --- 绘图与保存 ---
    # Figure 1: Volumetric stress 相位调制图
    plt.figure(1, figsize=(60, 16))
    RPM_raVol, ph_shiftVol, _, _ = mf.modulation_phase(AllphaVol[:, 1] * 180 / np.pi, bin_phase,
                                                       AllphaVol_ref[:, 1] * 180 / np.pi, 'yes3')
    plt.savefig(config.TM_all_region_Vol, format='pdf')

    # Figure 2: Normal stress 相位调制图
    plt.figure(2, figsize=(60, 16))
    RPM_raN, ph_shiftN, _, _ = mf.modulation_phase(AllphaN[:, 1] * 180 / np.pi, bin_phase,
                                                     AllphaN_ref[:, 1] * 180 / np.pi, 'yes3')
    plt.savefig(config.TM_all_region_N, format='pdf')

    # Figure 3: Shear stress 相位调制图
    plt.figure(3, figsize=(60, 16))
    RPM_raS, ph_shiftS, _, _ = mf.modulation_phase(AllphaS[:, 1] * 180 / np.pi, bin_phase,
                                                     AllphaS_ref[:, 1] * 180 / np.pi, 'yes3')
    plt.savefig(config.TM_all_region_S, format='pdf')

    # Figure 4: Coulomb stress 相位调制图
    plt.figure(4, figsize=(60, 16))
    RPM_raCFS, ph_shiftCFS, _, _ = mf.modulation_phase(AllphaCFS[:, 1] * 180 / np.pi, bin_phase,
                                                       AllphaCFS_ref[:, 1] * 180 / np.pi, 'yes3')
    plt.savefig(config.TM_all_region_CFS, format='pdf')



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

    fig =vis.plot_tidal_sensitivity_2x2(
        # Volumetric (第一个子图)
        a_estimated_Vol, c_estimated_Vol, delta_a_Vol, delta_c_Vol,
        Stress_AM_bk_Vol, Stress_AM_Vol, shear_stress_Vol, shear_stress_kPa_Vol, event_rate_Vol,
        # Normal (第二个子图)
        a_estimated_N, c_estimated_N, delta_a_N, delta_c_N,
        Stress_AM_bk_N, Stress_AM_N, shear_stress_N, shear_stress_kPa_N, event_rate_N,
        # Shear (第三个子图)
        a_estimated_S, c_estimated_S, delta_a_S, delta_c_S,
        Stress_AM_bk_S, Stress_AM_S, shear_stress_S, shear_stress_kPa_S, event_rate_S,
        # Coulomb (第四个子图)
        a_estimated_CFS, c_estimated_CFS, delta_a_CFS, delta_c_CFS,
        Stress_AM_bk_CFS, Stress_AM_CFS, shear_stress_CFS, shear_stress_kPa_CFS, event_rate_CFS
    )
    # 保存图片到配置文件中指定的路径
    fig.savefig(config.TM_all_region_tidal_sensitivity_cfs, format='pdf')
    plt.close(fig)
