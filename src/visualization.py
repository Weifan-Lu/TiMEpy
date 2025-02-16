import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from matplotlib import rcParams
# Set global font to Arial
rcParams['font.family'] = 'Arial'


def plot_p_value_subplot(subplot_idx, y_data, time_series, komo_t, title_text, ax=None):
    """
    绘制 p-value 子图（对数 y 轴）
    """
    if ax is None:
        ax = plt.gca()
    ax.plot(time_series, y_data, color='royalblue', marker='o', markersize=6,
            markerfacecolor='royalblue', markeredgecolor='black', linewidth=1.5, alpha=0.8,linestyle='-')
    ax.axhline(0.05, color='k', linestyle='--', linewidth=1.5)
    ax.axvline(komo_t, color='r', linewidth=1.5)
    ax = set_common_elements(ax, time_series, 'Time', 'P-value', title_text)
    # 设置 y 轴为对数刻度
    ax.set_yscale('log')
    return ax


def plot_amp_value_subplot(subplot_idx, y_data, time_series, komo_t, title_text, ax=None):
    """
    绘制振幅子图
    """
    if ax is None:
        ax = plt.gca()
    ax.plot(time_series, y_data, color='royalblue', marker='o', markersize=6,
            markerfacecolor='royalblue', markeredgecolor='black', linewidth=1.5, alpha=0.8,linestyle='-')
    ax.axvline(komo_t, color='r', linewidth=1.5)
    ax = set_common_elements(ax, time_series, 'Time', 'Amplitude', title_text)
    return ax


def plot_pha_value_subplot(subplot_idx, y_data, time_series, komo_t, title_text, ax=None):
    """
    绘制相位偏移子图
    """
    if ax is None:
        ax = plt.gca()
    ax.plot(time_series, y_data, color='royalblue', marker='o', markersize=6,
            markerfacecolor='royalblue', markeredgecolor='black', linewidth=1.5, alpha=0.8,linestyle='-')
    ax.axvline(komo_t, color='r', linewidth=1.5)
    ax = set_common_elements(ax, time_series, 'Time', 'Phase (°)', title_text)
    return ax


def plot_a_value_subplot(subplot_idx, A_est_CFS, time_series, komo_t, title_text, ax=None):
    """
    绘制灵敏度子图，假定 A_est_CFS 为 numpy 数组，其中第 1 列为数值，第 3 列为误差
    """
    if ax is None:
        ax = plt.gca()
    A_est_CFS = np.array(A_est_CFS)
    # 优化 errorbar
    ax.errorbar(
        time_series, A_est_CFS[:, 0], yerr=A_est_CFS[:, 2], fmt='o',
        markersize=6, markerfacecolor='royalblue', markeredgecolor='black',  # 美化标记
        capsize=4, capthick=1.2, elinewidth=1, alpha=0.8, color='royalblue',  # 美化误差棒
        label='S',linestyle='-')

    ax.axvline(komo_t, color='r', linewidth=1.5)
    ax = set_common_elements(ax, time_series, 'Time', 'Tidal sensitivity (α)', title_text)
    return ax


def set_common_elements(ax, time_series, x_label, y_label, title_text):
    xlim_range = [min(time_series), max(time_series) + 365 * 2]
    start_year = mdates.num2date(xlim_range[0]).year
    end_year = mdates.num2date(xlim_range[1]).year
    years = np.arange(start_year, end_year + 1, 2)
    year_ticks = [mdates.date2num(datetime(year, 1, 1)) for year in years]
    ax.set_xticks(year_ticks)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(title_text, fontsize=16)
    ax.grid(True)
    ax.tick_params(axis='both', labelsize=12)
    return ax


def plot_phase_modulation(ph1, Prob_o, Po, G, p, wbin, PM_ra, ph_shift):
    """
    该函数基于给定的数据生成三个子图
    """
    bar_color = [0.8, 0.8, 0.8]
    bar_edge_color = 'k'
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))

    # 第一个子图：观测分布
    axs[0].bar(ph1, Prob_o, width=wbin, color=bar_color, edgecolor=bar_edge_color)
    set_graph_elements(axs[0], 'Phase (degree)', 'Probability Density', 'P_obs (Observed distribution)')

    # 第二个子图：参考分布
    axs[1].bar(ph1, Po, width=wbin, color=bar_color, edgecolor=bar_edge_color)
    set_graph_elements(axs[1], 'Phase (degree)', 'Probability Density', 'P_ref (Uniform distribution)')

    # 第三个子图：归一化比值与拟合曲线
    axs[2].bar(ph1, Prob_o / Po, width=wbin, color=bar_color, edgecolor=bar_edge_color)
    axs[2].plot(ph1, np.ones_like(ph1), 'k--', linewidth=2)  # 参考线 y=1
    axs[2].plot(ph1, (G @ p + 1), 'r-', linewidth=2)
    axs[2].set_xlabel('Phase (degree)')
    axs[2].set_ylabel('Normalized Density')
    axs[2].set_title(f'P_obs/P_ref: {PM_ra:.2f}, Phase shift: {ph_shift:.1f}°')
    axs[2].set_xlim([-180, 180])
    axs[2].set_xticks(np.arange(-180, 181, 45))
    axs[2].grid(True)

    plt.tight_layout()


def set_graph_elements(ax, xlabel, ylabel, title):
    """
    设置图形元素
    """
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim([-180, 180])
    ax.set_xticks(np.arange(-180, 181, 45))
    ax.grid(True)


def plot_modulation_stress(Stress_AM_bk, Stress_AM, shear_stress, shear_stress_kPa, event_rate, a_estimated,
                           C_estimated, delta_a, delta_c):
    """
    绘制应力数据直方图以及拟合结果图。

    参数:
        Stress_AM_bk: 参考应力数据（一维数组）
        Stress_AM: 观测应力数据（一维数组）
        shear_stress: 分箱后应力的中心值（经过非零过滤）
        shear_stress_kPa: shear_stress 转换为 kPa 后的值
        event_rate: 每个箱中的事件比率
        a_estimated: 优化得到的 a 参数
        C_estimated: 优化得到的 C 参数
        delta_a: a 的 95% 置信区间估计
        delta_c: C 的 95% 置信区间估计

    返回:
        fig, ax1, ax2: 绘图对象，方便在外部进行显示或保存
    """
    fig, ax1 = plt.subplots()

    # 绘制参考应力数据和观测应力数据的直方图
    counts1, edges1 = np.histogram(Stress_AM_bk, bins=40)
    counts1 = counts1 / np.max(counts1)
    centers1 = (edges1[:-1] + edges1[1:]) / 2
    counts2, _ = np.histogram(Stress_AM, bins=edges1)
    counts2 = counts2 / np.max(counts2)

    ax1.bar(centers1, counts1, width=np.diff(edges1), color='none', edgecolor='k',
            linewidth=1.5, label='Stress AM bk')
    ax1.bar(centers1, counts2, width=np.diff(edges1), color='none', edgecolor='r',
            linewidth=1.5, label='Stress AM')

    ax1.set_xlabel('Shear Stress (Pa)')
    ax1.set_ylabel('Normalized frequency')

    # 绘制观测数据点和拟合曲线
    ax2 = ax1.twinx()
    fitted_rate = C_estimated * np.exp(a_estimated * shear_stress_kPa)
    ax2.plot(shear_stress, event_rate, 'ko', linewidth=1.5)
    ax2.plot(shear_stress, fitted_rate, 'r-', linewidth=2, label='Fitted')

    ax2.set_ylabel('Tremor Rate (events/hour)')
    ax1.tick_params(axis='both', labelsize=16)
    ax1.legend()

    # 添加参数估计值及误差
    x_pos = min(Stress_AM_bk) * 0.9  # 确保 x 轴范围合适
    ax1.text(x_pos, 0.95, r'$\alpha = {:.2f} \pm {:.3f} \,\mathrm{{kPa}}^{{-1}}$'.format(a_estimated, delta_a),
             fontsize=14, color='blue')
    ax1.text(x_pos, 0.9, r'$C = {:.4f} \pm {:.3f} \,\mathrm{{h}}^{{-1}}$'.format(C_estimated, delta_c),
             fontsize=14, color='blue')

    # 注意：这里不调用 plt.show()，由外部统一管理显示或保存
    return fig, ax1, ax2

def plot_tidal_sensitivity_2x2(
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
):
    """
    绘制 2x2 子图：
        第1个子图：Volumetric Strain
        第2个子图：Normal Stress
        第3个子图：Shear Stress
        第4个子图：Coulomb Stress

    每个子图中：
      - 左侧 y 轴绘制参考应力数据（Stress_AM_bk）与观测应力数据（Stress_AM）的归一化直方图；
      - 右侧 y 轴绘制拟合数据点（event_rate）及拟合曲线（fitted_rate = C * exp(a * shear_stress_kPa)）；
      - 在左上角标注参数值（α 与 C 及其 95% 置信区间）。

    绘图结果保存为 PDF 文件，保存路径由 config.TM_all_region_tidal_sensitivity_cfs 指定。
    """

    # 创建 2x2 子图
    print(len(shear_stress_Vol), len(event_rate_Vol))
    print(len(shear_stress_S), len(event_rate_S))
    print(len(shear_stress_N), len(event_rate_N))
    print(len(shear_stress_CFS), len(event_rate_CFS))

    fig, axs = plt.subplots(2, 2, figsize=(20, 16))
    plt.rcParams.update({'font.size': 24, 'font.family': 'arial'})  # 设置字体大小和样式
    n_front = 24
    n_front_lab = 20
    n_linewidth = 2
    ###############################################
    # 子图 1：Volumetric Strain
    ###############################################
    ax = axs[0, 0]
    counts1, edges1 = np.histogram(Stress_AM_bk_Vol, bins=20)
    counts1 = counts1 / np.max(counts1)
    centers1 = (edges1[:-1] + edges1[1:]) / 2
    counts2, _ = np.histogram(Stress_AM_Vol, bins=edges1)
    counts2 = counts2 / np.max(counts2)

    ax.bar(centers1, counts1, width=np.diff(edges1), color='none', edgecolor='k',
           linewidth=n_linewidth, label='Stress AM bk', alpha=0.7)
    ax.bar(centers1, counts2, width=np.diff(edges1), color='none', edgecolor='r',
           linewidth=n_linewidth, label='Stress AM', alpha=0.7)
    ax.set_xlabel(r'Volumetric Strain ($\times 10^{10}$)', fontsize=n_front,fontname='Arial')
    ax.set_ylabel('Normalized frequency', fontsize=n_front)
    # ax.legend(fontsize=14)
    ax.tick_params(axis='both', labelsize=n_front_lab)
    ax.grid(True, linestyle='--', alpha=0.6)

    ax2 = ax.twinx()
    fitted_rate_Vol = c_estimated_Vol * np.exp(a_estimated_Vol * shear_stress_kPa_Vol)
    ax2.plot(shear_stress_Vol, event_rate_Vol, 'ko', markersize=10, label='Data',mfc='none',mew=2)
    ax2.plot(shear_stress_Vol, fitted_rate_Vol, 'r-', linewidth=2.5, label='Fitted')
    ax2.set_ylabel('Event rate (events/hour)', fontsize=n_front)
    ax2.tick_params(axis='both', labelsize=n_front_lab)

    # 使用归一化坐标添加文本标注
    ax.text(0.05, 0.95, r'$\alpha = {:.2f} \pm {:.3f}\,\mathrm{{kPa}}^{{-1}}$'
            .format(a_estimated_Vol, delta_a_Vol),
            fontsize=16, color='blue', transform=ax.transAxes)
    ax.text(0.05, 0.90, r'$C = {:.4f} \pm {:.3f}\,\mathrm{{h}}^{{-1}}$'
            .format(c_estimated_Vol, delta_c_Vol),
            fontsize=16, color='blue', transform=ax.transAxes)

    ###############################################
    # 子图 2：Normal Stress
    ###############################################
    ax = axs[0, 1]
    counts1, edges1 = np.histogram(Stress_AM_bk_N, bins=20)
    counts1 = counts1 / np.max(counts1)
    centers1 = (edges1[:-1] + edges1[1:]) / 2
    counts2, _ = np.histogram(Stress_AM_N, bins=edges1)
    counts2 = counts2 / np.max(counts2)

    ax.bar(centers1, counts1, width=np.diff(edges1), color='none', edgecolor='k',
           linewidth=n_linewidth, label='Stress AM bk', alpha=0.7)
    ax.bar(centers1, counts2, width=np.diff(edges1), color='none', edgecolor='r',
           linewidth=n_linewidth, label='Stress AM', alpha=0.7)
    ax.set_xlabel('Normal stress (Pa)', fontsize=n_front)
    ax.set_ylabel('Normalized frequency', fontsize=n_front)
    # ax.legend(fontsize=14)
    ax.tick_params(axis='both', labelsize=n_front_lab)
    ax.grid(True, linestyle='--', alpha=0.6)

    ax2 = ax.twinx()
    fitted_rate_N = c_estimated_N * np.exp(a_estimated_N * shear_stress_kPa_N)
    ax2.plot(shear_stress_N, event_rate_N, 'ko', markersize=10, label='Data',mfc='none',mew=2)
    ax2.plot(shear_stress_N, fitted_rate_N, 'r-', linewidth=2.5, label='Fitted')
    ax2.set_ylabel('Event rate (events/hour)', fontsize=n_front)
    ax2.tick_params(axis='both', labelsize=n_front_lab)

    ax.text(0.05, 0.95, r'$\alpha = {:.2f} \pm {:.3f}\,\mathrm{{kPa}}^{{-1}}$'
            .format(a_estimated_N, delta_a_N),
            fontsize=16, color='blue', transform=ax.transAxes)
    ax.text(0.05, 0.90, r'$C = {:.4f} \pm {:.3f}\,\mathrm{{h}}^{{-1}}$'
            .format(c_estimated_N, delta_c_N),
            fontsize=16, color='blue', transform=ax.transAxes)

    ###############################################
    # 子图 3：Shear Stress
    ###############################################
    ax = axs[1, 0]
    counts1, edges1 = np.histogram(Stress_AM_bk_S, bins=20)
    counts1 = counts1 / np.max(counts1)
    centers1 = (edges1[:-1] + edges1[1:]) / 2
    counts2, _ = np.histogram(Stress_AM_S, bins=edges1)
    counts2 = counts2 / np.max(counts2)

    ax.bar(centers1, counts1, width=np.diff(edges1), color='none', edgecolor='k',
           linewidth=n_linewidth, label='Stress AM bk', alpha=0.7)
    ax.bar(centers1, counts2, width=np.diff(edges1), color='none', edgecolor='r',
           linewidth=n_linewidth, label='Stress AM', alpha=0.7)
    ax.set_xlabel('Shear stress (Pa)', fontsize=n_front)
    ax.set_ylabel('Normalized frequency', fontsize=n_front)
    # ax.legend(fontsize=14)
    ax.tick_params(axis='both', labelsize=n_front_lab)
    ax.grid(True, linestyle='--', alpha=0.6)

    ax2 = ax.twinx()
    fitted_rate_S = c_estimated_S * np.exp(a_estimated_S * shear_stress_kPa_S)
    ax2.plot(shear_stress_S, event_rate_S, 'ko', markersize=10, label='Data',mfc='none',mew=2)
    ax2.plot(shear_stress_S, fitted_rate_S, 'r-', linewidth=2.5, label='Fitted')
    ax2.set_ylabel('Event rate (events/hour)', fontsize=n_front)
    ax2.tick_params(axis='both', labelsize=n_front_lab)

    ax.text(0.05, 0.95, r'$\alpha = {:.2f} \pm {:.3f}\,\mathrm{{kPa}}^{{-1}}$'
            .format(a_estimated_S, delta_a_S),
            fontsize=16, color='blue', transform=ax.transAxes)
    ax.text(0.05, 0.90, r'$C = {:.4f} \pm {:.3f}\,\mathrm{{h}}^{{-1}}$'
            .format(c_estimated_S, delta_c_S),
            fontsize=16, color='blue', transform=ax.transAxes)

    ###############################################
    # 子图 4：Coulomb Stress
    ###############################################
    ax = axs[1, 1]
    counts1, edges1 = np.histogram(Stress_AM_bk_CFS, bins=20)
    counts1 = counts1 / np.max(counts1)
    centers1 = (edges1[:-1] + edges1[1:]) / 2
    counts2, _ = np.histogram(Stress_AM_CFS, bins=edges1)
    counts2 = counts2 / np.max(counts2)

    ax.bar(centers1, counts1, width=np.diff(edges1), color='none', edgecolor='k',
           linewidth=n_linewidth, label='Stress AM bk', alpha=0.7)
    ax.bar(centers1, counts2, width=np.diff(edges1), color='none', edgecolor='r',
           linewidth=n_linewidth, label='Stress AM', alpha=0.7)
    ax.set_xlabel('Coulomb stress (Pa)', fontsize=n_front)
    ax.set_ylabel('Normalized frequency', fontsize=n_front)
    # ax.legend(fontsize=14)
    ax.tick_params(axis='both', labelsize=n_front_lab)
    ax.grid(True, linestyle='--', alpha=0.6)

    ax2 = ax.twinx()
    fitted_rate_CFS = c_estimated_CFS * np.exp(a_estimated_CFS * shear_stress_kPa_CFS)
    ax2.plot(shear_stress_CFS, event_rate_CFS, 'ko', markersize=10, label='Data',mfc='none',mew=2)
    ax2.plot(shear_stress_CFS, fitted_rate_CFS, 'r-', linewidth=2.5, label='Fitted')
    ax2.set_ylabel('Event rate (events/hour)', fontsize=16)
    ax2.tick_params(axis='both', labelsize=n_front_lab)

    ax.text(0.05, 0.95, r'$\alpha = {:.2f} \pm {:.3f}\,\mathrm{{kPa}}^{{-1}}$'
            .format(a_estimated_CFS, delta_a_CFS),
            fontsize=16, color='blue', transform=ax.transAxes)
    ax.text(0.05, 0.90, r'$C = {:.4f} \pm {:.3f}\,\mathrm{{h}}^{{-1}}$'
            .format(c_estimated_CFS, delta_c_CFS),
            fontsize=16, color='blue', transform=ax.transAxes)

    ###############################################
    # 调整布局并保存图形
    ###############################################
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.3, hspace=0.3)  # 调整子图间距
    return fig

def plot_stress_earthquake(t_decluster, volumetric_strain, normal_stress, shear_stress, stress_cfs,
                           t_tide, StressVol, StressN, StressS, StressCFS, config):
    # 绘图
    t_start_dt = config.start_time
    t_end_dt = config.end_time
    fig, axes = plt.subplots(4, 1, figsize=(12, 12))

    # Volumetric strain
    axes[0].plot(t_tide, StressVol, color='black', linestyle='-', linewidth=0.5, label='Volumetric strain')
    axes[0].plot(t_decluster, volumetric_strain, 'ro', markersize=2)
    axes[0].set_ylabel('Volumetric strain', fontsize=16)
    axes[0].grid()

    # Normal stress
    axes[1].plot(t_tide, StressN, color='black', linestyle='-', linewidth=0.5, label='Normal stress')
    axes[1].plot(t_decluster, normal_stress, 'ro', markersize=2)
    axes[1].set_ylabel('Normal stress (Pa)', fontsize=16)
    axes[1].grid()

    # Shear stress
    axes[2].plot(t_tide, StressS, color='black', linestyle='-', linewidth=0.5, label='Shear stress')
    axes[2].plot(t_decluster, shear_stress, 'ro', markersize=2)
    axes[2].set_ylabel('Shear stress (Pa)', fontsize=16)
    axes[2].grid()

    # Coulomb stress
    axes[3].plot(t_tide, StressCFS, color='black', linestyle='-', linewidth=0.5, label='Coulomb stress')
    axes[3].plot(t_decluster, stress_cfs, 'ro', markersize=2)
    axes[3].set_ylabel('Coulomb stress (Pa)', fontsize=16)
    axes[3].grid()

    # 格式化 x 轴：只在最后一个子图显示 x 轴标签和刻度
    for i, ax in enumerate(axes):
        ax.set_xlim([t_start_dt, t_end_dt])
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        ax.tick_params(axis='y', labelsize=16)
        if i < len(axes) - 1:
            ax.tick_params(axis='x', which='both', labelbottom=False)  # 隐藏 x 轴标签
        else:
            ax.tick_params(axis='x', rotation=45, labelsize=16)
            ax.set_xlabel('Time', fontsize=16)

    plt.tight_layout()
    plt.savefig(config.TM_all_region_stress, format='pdf', dpi=300)
    # plt.show()
