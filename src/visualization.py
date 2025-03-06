import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from matplotlib import rcParams, cm, colors
from matplotlib.colors import ListedColormap
from matplotlib.cm import ScalarMappable

# rcParams['font.family'] = 'Arial'
# rcParams['font.family'] = 'DejaVu Sans'

def plot_p_value_subplot(subplot_idx, y_data, time_series, komo_t, title_text, ax=None, font_size=14):
    """
    Plot p-value subplot (logarithmic y-axis)
    """
    if ax is None:
        ax = plt.gca()
    ax.plot(time_series, y_data, color='royalblue', marker='o', markersize=6,
            markerfacecolor='royalblue', markeredgecolor='black', linewidth=2,
            alpha=0.8, linestyle='-')
    ax.axhline(0.05, color='k', linestyle='--', linewidth=2)
    ax.axvline(komo_t, color='r', linewidth=2)
    # Assuming set_common_elements can accept font_size parameter to set title, x, and y labels
    ax = set_common_elements(ax, time_series, 'Time', 'P-value', title_text, font_size=font_size)
    ax.set_yscale('log')
    return ax

def plot_amp_value_subplot(subplot_idx, y_data, time_series, komo_t, title_text, ax=None, font_size=14):
    """
    Plot amplitude subplot
    """
    if ax is None:
        ax = plt.gca()
    ax.plot(time_series, y_data, color='royalblue', marker='o', markersize=6,
            markerfacecolor='royalblue', markeredgecolor='black', linewidth=2,
            alpha=0.8, linestyle='-')
    ax.axvline(komo_t, color='r', linewidth=2)
    ax = set_common_elements(ax, time_series, 'Time', 'Phase amplitude ($P$)', title_text, font_size=font_size)
    return ax

def plot_pha_value_subplot(subplot_idx, y_data, time_series, komo_t, title_text, ax=None, font_size=14):
    """
    Plot phase shift subplot
    """
    if ax is None:
        ax = plt.gca()
    ax.plot(time_series, y_data, color='royalblue', marker='o', markersize=6,
            markerfacecolor='royalblue', markeredgecolor='black', linewidth=2,
            alpha=0.8, linestyle='-')
    ax.axvline(komo_t, color='r', linewidth=2)
    # ax.set_ylim([-180,180])
    ax = set_common_elements(ax, time_series, 'Time', 'Phase shift (°)', title_text, font_size=font_size)
    return ax

def plot_a_value_subplot(subplot_idx, A_est_CFS, time_series, komo_t, title_text, ax=None, font_size=14):
    """
    Plot sensitivity subplot, assuming A_est_CFS is a numpy array where the first column contains values
    and the third column contains errors
    """
    if ax is None:
        ax = plt.gca()
    A_est_CFS = np.array(A_est_CFS)
    ax.errorbar(
        time_series, A_est_CFS[:, 0], yerr=A_est_CFS[:, 2], fmt='o',
        markersize=6, markerfacecolor='royalblue', markeredgecolor='black',
        capsize=4, capthick=1.2, elinewidth=2, alpha=0.8, color='royalblue',
        label='S', linestyle='-')
    ax.axvline(komo_t, color='r', linewidth=2)
    ax = set_common_elements(ax, time_series, 'Time', 'Tidal sensitivity (α)', title_text, font_size=font_size)
    return ax

def set_common_elements(ax, time_series, x_label, y_label, title_text, font_size):
    xlim_range = [min(time_series), max(time_series) + 365 * 1.5]
    start_year = mdates.num2date(xlim_range[0]).year
    end_year = mdates.num2date(xlim_range[1]).year
    years = np.arange(start_year, end_year + 1, 2)
    year_ticks = [mdates.date2num(datetime(year, 1, 1)) for year in years]
    ax.set_xticks(year_ticks)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

    ax.set_xlabel(x_label, fontsize=font_size)
    ax.set_ylabel(y_label, fontsize=font_size)
    ax.set_title(title_text, fontsize=font_size)
    ax.grid(True)
    ax.tick_params(axis='both', labelsize=font_size)
    return ax

def plot_phase_modulation(ph1, Prob_o, Po, wbin, PM_ra, ph_shift, model, phi_dense):
    """
    This function generates three subplots based on the given data.

    Parameters:
        ph1: X-axis data (phase in degrees)
        Prob_o: Observed event count (or frequency)
        Po: Reference event count (or frequency)
        wbin: Bin width for the bars
        PM_ra: Numeric value displayed in the title (rounded to two decimal places)
        ph_shift: Phase shift displayed in the title (rounded to one decimal place)
        model: Function for generating the fitted curve
        phi_dense: X-axis data for the fitted curve

    Returns:
        fig, axs: Matplotlib figure and subplot objects
    """
    n_front = 16
    bar_color = [0.8, 0.8, 0.8]
    bar_edge_color = 'k'
    fig, axs = plt.subplots(1, 3, figsize=(24,6))

    # First subplot: Observed distribution
    axs[0].bar(ph1, Prob_o, width=wbin, color=bar_color, edgecolor=bar_edge_color)
    axs[0].set_xlabel(r'Phase ($^\circ$)', fontsize=n_front)
    axs[0].set_ylabel('Probability Density', fontsize=n_front)
    axs[0].set_title('P_obs (Observed distribution)', fontsize=n_front)
    axs[0].tick_params(axis='both', labelsize=n_front)

    # Second subplot: Reference distribution
    axs[1].bar(ph1, Po, width=wbin, color=bar_color, edgecolor=bar_edge_color)
    axs[1].set_xlabel(r'Phase ($^\circ$)', fontsize=n_front)
    axs[1].set_ylabel('Probability Density', fontsize=n_front)
    axs[1].set_title('P_ref (Uniform distribution)', fontsize=n_front)
    axs[1].tick_params(axis='both', labelsize=n_front)

    # Third subplot: Normalized ratio and fitted curve
    axs[2].bar(ph1, Prob_o / Po, width=wbin, color=bar_color, edgecolor=bar_edge_color)
    axs[2].plot(ph1, np.ones_like(ph1), 'k--', linewidth=2)  # Reference line y=1
    axs[2].plot(phi_dense, model(phi_dense) + 1, 'r-', label='Fitted', linewidth=2)
    axs[2].set_xlabel(r'Phase ($^\circ$)', fontsize=n_front)
    axs[2].set_ylabel('$P_{obs}/P_{ref}$', fontsize=n_front)
    axs[2].set_title(f'Phase amplitude: {PM_ra:.2f}, Phase shift: {ph_shift:.1f}°', fontsize=n_front)
    axs[2].set_xlim([-180, 180])
    ticks = np.arange(-180, 181, 45)
    axs[2].set_xticks(ticks)
    axs[2].set_xticklabels(ticks, fontsize=n_front)
    axs[2].tick_params(axis='both', labelsize=n_front)

def set_graph_elements(ax, xlabel, ylabel, title):
    """
    Set graph elements such as labels, title, and grid.
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
    Plot histograms of stress data and fitting results.

    Parameters:
        Stress_AM_bk: Reference stress data (1D array)
        Stress_AM: Observed stress data (1D array)
        shear_stress: Center values of binned stress (filtered for nonzero values)
        shear_stress_kPa: shear_stress converted to kPa
        event_rate: Event rate in each bin
        a_estimated: Optimized parameter a
        C_estimated: Optimized parameter C
        delta_a: 95% confidence interval estimate for a
        delta_c: 95% confidence interval estimate for C

    Returns:
        fig, ax1, ax2: Plot objects for external display or saving
    """
    # Set a uniform font size
    font_size = 12

    fig, ax1 = plt.subplots(figsize=(6, 5))

    # Plot histograms of reference stress data and observed stress data
    edges1 = np.linspace(np.min(Stress_AM), np.max(Stress_AM), 20)
    counts1, _ = np.histogram(Stress_AM_bk, bins=edges1)
    counts1 = counts1 / np.sum(counts1)
    centers1 = (edges1[:-1] + edges1[1:]) / 2
    counts2, _ = np.histogram(Stress_AM, bins=edges1)
    counts2 = counts2 / np.sum(counts2)

    ax1.bar(centers1, counts1, width=np.diff(edges1), color=[0.8, 0.8, 0.8], edgecolor='k',
            linewidth=1.5, label='Ref', alpha=0.7)
    ax1.bar(centers1, counts2, width=np.diff(edges1), color='none', edgecolor='r',
            linewidth=1.5, label='Obs')

    ax1.set_xlabel('Coulomb stress (Pa)', fontsize=font_size)
    ax1.set_ylabel('Frequency', fontsize=font_size)

    # Plot observed data points and fitted curve
    ax2 = ax1.twinx()
    fitted_rate = C_estimated * np.exp(a_estimated * shear_stress_kPa)
    ax2.plot(shear_stress, event_rate, 'ko', markersize=6, label='Data', mfc='none', mew=1)
    ax2.plot(shear_stress, fitted_rate, 'r-', linewidth=2, label='Fitted')
    ax2.set_ylabel('Obs. / Ref. rate', fontsize=font_size)

    # Set uniform tick label size
    ax1.tick_params(axis='both', labelsize=font_size)
    ax2.tick_params(axis='both', labelsize=font_size)

    # Add legend (position adjustable as needed)
    leg1 = ax1.legend(fontsize=font_size, loc='upper left')

    x_pos = max(Stress_AM_bk) * 0.3
    y_pos1 = np.max(event_rate) * 0.95  # Adjust based on data
    y_pos2 = np.max(event_rate) * 0.90

    ax2.text(x_pos, y_pos1,
             r'$\alpha = {:.2f} \pm {:.3f}\,\mathrm{{kPa}}^{{-1}}$'.format(a_estimated, delta_a),
             fontsize=font_size, color='blue')
    ax2.text(x_pos, y_pos2,
             r'$C = {:.4f} \pm {:.3f}\,\mathrm{{h}}^{{-1}}$'.format(C_estimated, delta_c),
             fontsize=font_size, color='blue')

    # plt.show()

    # Note: plt.show() is not called here, allowing external control for display or saving
    return fig, ax1, ax2


def plot_tidal_sensitivity_2x2(
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
):
    """
        Plot 2x2 subplots:
            1st subplot: Volumetric Strain
            2nd subplot: Normal Stress
            3rd subplot: Shear Stress
            4th subplot: Coulomb Stress

        In each subplot:
          - The left y-axis shows the normalized histograms of reference stress data (Stress_AM_bk)
            and observed stress data (Stress_AM).
          - The right y-axis shows the fitted data points (event_rate) and the fitted curve
            (fitted_rate = C * exp(a * shear_stress_kPa)).
          - The top-left corner displays parameter values (α and C with their 95% confidence intervals).

        The plot is saved as a PDF file, with the save path specified by config.TM_all_region_tidal_sensitivity_cfs.
        """

    # 创建 2x2 子图
    # print(len(shear_stress_Vol), len(event_rate_Vol))
    # print(len(shear_stress_S), len(event_rate_S))
    # print(len(shear_stress_N), len(event_rate_N))
    # print(len(shear_stress_CFS), len(event_rate_CFS))

    fig, axs = plt.subplots(2, 2, figsize=(20, 16))
    # plt.rcParams.update({'font.size': 24, 'font.family': 'arial'})
    n_front = 24
    n_front_lab = 20
    n_linewidth = 2
    ###############################################
    # 1：Volumetric Strain
    ###############################################
    ax = axs[0, 0]
    edges1 = np.linspace(np.min(Stress_AM_Vol),np.max(Stress_AM_Vol), 20)
    counts1, _ = np.histogram(Stress_AM_bk_Vol, bins=edges1)
    counts1 = counts1 / np.sum(counts1)
    centers1 = (edges1[:-1] + edges1[1:]) / 2
    counts2, _ = np.histogram(Stress_AM_Vol, bins=edges1)
    counts2 = counts2 / np.sum(counts2)

    ax.bar(centers1, counts1, width=np.diff(edges1), color=[0.8, 0.8, 0.8], edgecolor='k',
           linewidth=n_linewidth, label='Ref', alpha=0.7)
    ax.bar(centers1, counts2, width=np.diff(edges1), color='none', edgecolor='r',
           linewidth=n_linewidth, label='Obs', alpha=0.7)
    ax.set_xlabel(r'Volumetric Strain ($\times 10^{10}$)', fontsize=n_front)
    ax.set_ylabel('Frequency', fontsize=n_front)
    ax.legend(fontsize=14)
    ax.tick_params(axis='both', labelsize=n_front_lab)
    ax.grid(True, linestyle='--', alpha=0.6)

    ax2 = ax.twinx()
    fitted_rate_Vol = c_estimated_Vol * np.exp(a_estimated_Vol * shear_stress_kPa_Vol)
    ax2.plot(shear_stress_Vol, event_rate_Vol, 'ko', markersize=10, label='Data',mfc='none',mew=2)
    ax2.plot(shear_stress_Vol, fitted_rate_Vol, 'r-', linewidth=2.5, label='Fitted')
    ax2.set_ylabel('Obs. / Ref. rate', fontsize=n_front)
    ax2.tick_params(axis='both', labelsize=n_front_lab)

    ax.text(0.05, 0.95, r'$\alpha = {:.2f} \pm {:.3f}\,\mathrm{{kPa}}^{{-1}}$'
            .format(a_estimated_Vol, delta_a_Vol),
            fontsize=n_front_lab, color='blue', transform=ax.transAxes)
    ax.text(0.05, 0.90, r'$C = {:.4f} \pm {:.3f}\,\mathrm{{h}}^{{-1}}$'
            .format(c_estimated_Vol, delta_c_Vol),
            fontsize=n_front_lab, color='blue', transform=ax.transAxes)

    ###############################################
    # 2：Normal Stress
    ###############################################
    ax = axs[0, 1]
    edges1 = np.linspace(np.min(Stress_AM_N),np.max(Stress_AM_N), 20)
    counts1, _ = np.histogram(Stress_AM_bk_N, bins=edges1)
    counts1 = counts1 / np.sum(counts1)
    centers1 = (edges1[:-1] + edges1[1:]) / 2
    counts2, _ = np.histogram(Stress_AM_N, bins=edges1)
    counts2 = counts2 / np.sum(counts2)

    ax.bar(centers1, counts1, width=np.diff(edges1), color=[0.8, 0.8, 0.8], edgecolor='k',
           linewidth=n_linewidth, label='Stress AM bk', alpha=0.7)
    ax.bar(centers1, counts2, width=np.diff(edges1), color='none', edgecolor='r',
           linewidth=n_linewidth, label='Stress AM', alpha=0.7)
    ax.set_xlabel('Normal stress (Pa)', fontsize=n_front)
    ax.set_ylabel('Frequency', fontsize=n_front)
    ax.tick_params(axis='both', labelsize=n_front_lab)
    ax.grid(True, linestyle='--', alpha=0.6)

    ax2 = ax.twinx()
    fitted_rate_N = c_estimated_N * np.exp(a_estimated_N * shear_stress_kPa_N)
    ax2.plot(shear_stress_N, event_rate_N, 'ko', markersize=10, label='Data',mfc='none',mew=2)
    ax2.plot(shear_stress_N, fitted_rate_N, 'r-', linewidth=2.5, label='Fitted')
    ax2.set_ylabel('Obs. / Ref. rate', fontsize=n_front)
    ax2.tick_params(axis='both', labelsize=n_front_lab)

    ax.text(0.05, 0.95, r'$\alpha = {:.2f} \pm {:.3f}\,\mathrm{{kPa}}^{{-1}}$'
            .format(a_estimated_N, delta_a_N),
            fontsize=n_front_lab, color='blue', transform=ax.transAxes)
    ax.text(0.05, 0.90, r'$C = {:.4f} \pm {:.3f}\,\mathrm{{h}}^{{-1}}$'
            .format(c_estimated_N, delta_c_N),
            fontsize=n_front_lab, color='blue', transform=ax.transAxes)

    ###############################################
    #3：Shear Stress
    ###############################################
    ax = axs[1, 0]
    edges1 = np.linspace(np.min(Stress_AM_S),np.max(Stress_AM_S), 20)
    counts1, _ = np.histogram(Stress_AM_bk_S, bins=edges1)
    counts1 = counts1 / np.sum(counts1)
    centers1 = (edges1[:-1] + edges1[1:]) / 2
    counts2, _ = np.histogram(Stress_AM_S, bins=edges1)
    counts2 = counts2 / np.sum(counts2)

    ax.bar(centers1, counts1, width=np.diff(edges1), color=[0.8, 0.8, 0.8], edgecolor='k',
           linewidth=n_linewidth, label='Stress AM bk', alpha=0.7)
    ax.bar(centers1, counts2, width=np.diff(edges1), color='none', edgecolor='r',
           linewidth=n_linewidth, label='Stress AM', alpha=0.7)
    ax.set_xlabel('Shear stress (Pa)', fontsize=n_front)
    ax.set_ylabel('Frequency', fontsize=n_front)
    ax.tick_params(axis='both', labelsize=n_front_lab)
    ax.grid(True, linestyle='--', alpha=0.6)

    ax2 = ax.twinx()
    fitted_rate_S = c_estimated_S * np.exp(a_estimated_S * shear_stress_kPa_S)
    ax2.plot(shear_stress_S, event_rate_S, 'ko', markersize=10, label='Data',mfc='none',mew=2)
    ax2.plot(shear_stress_S, fitted_rate_S, 'r-', linewidth=2.5, label='Fitted')
    ax2.set_ylabel('Obs. / Ref. rate', fontsize=n_front)
    ax2.tick_params(axis='both', labelsize=n_front_lab)

    ax.text(0.05, 0.95, r'$\alpha = {:.2f} \pm {:.3f}\,\mathrm{{kPa}}^{{-1}}$'
            .format(a_estimated_S, delta_a_S),
            fontsize=n_front_lab, color='blue', transform=ax.transAxes)
    ax.text(0.05, 0.90, r'$C = {:.4f} \pm {:.3f}\,\mathrm{{h}}^{{-1}}$'
            .format(c_estimated_S, delta_c_S),
            fontsize=n_front_lab, color='blue', transform=ax.transAxes)

    ###############################################
    #4：Coulomb Stress
    ###############################################
    ax = axs[1, 1]
    edges1 = np.linspace(np.min(Stress_AM_CFS),np.max(Stress_AM_CFS), 20)
    counts1, _ = np.histogram(Stress_AM_bk_CFS, bins=edges1)
    counts1 = counts1 / np.sum(counts1)
    centers1 = (edges1[:-1] + edges1[1:]) / 2
    counts2, _ = np.histogram(Stress_AM_CFS, bins=edges1)
    counts2 = counts2 / np.sum(counts2)

    ax.bar(centers1, counts1, width=np.diff(edges1), color=[0.8, 0.8, 0.8], edgecolor='k',
           linewidth=n_linewidth, label='Stress AM bk', alpha=0.7)
    ax.bar(centers1, counts2, width=np.diff(edges1), color='none', edgecolor='r',
           linewidth=n_linewidth, label='Stress AM', alpha=0.7)
    ax.set_xlabel('Coulomb stress (Pa)', fontsize=n_front)
    ax.set_ylabel('Frequency', fontsize=n_front)
    ax.tick_params(axis='both', labelsize=n_front_lab)
    ax.grid(True, linestyle='--', alpha=0.6)

    ax2 = ax.twinx()
    fitted_rate_CFS = c_estimated_CFS * np.exp(a_estimated_CFS * shear_stress_kPa_CFS)
    ax2.plot(shear_stress_CFS, event_rate_CFS, 'ko', markersize=10, label='Data',mfc='none',mew=2)
    ax2.plot(shear_stress_CFS, fitted_rate_CFS, 'r-', linewidth=2.5, label='Fitted')
    ax2.set_ylabel('Obs. / Ref. rate', fontsize=n_front)
    ax2.tick_params(axis='both', labelsize=n_front_lab)

    ax.text(0.05, 0.95, r'$\alpha = {:.2f} \pm {:.3f}\,\mathrm{{kPa}}^{{-1}}$'
            .format(a_estimated_CFS, delta_a_CFS),
            fontsize=n_front_lab, color='blue', transform=ax.transAxes)
    ax.text(0.05, 0.90, r'$C = {:.4f} \pm {:.3f}\,\mathrm{{h}}^{{-1}}$'
            .format(c_estimated_CFS, delta_c_CFS),
            fontsize=n_front_lab, color='blue', transform=ax.transAxes)

    ###############################################
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    return fig


def plot_stress_earthquake(t_decluster, volumetric_strain, normal_stress, shear_stress, stress_cfs,
                           t_tide, StressVol, StressN, StressS, StressCFS, config):
    # Plot earthquake stress data

    # Retrieve the start and end times from the configuration
    t_start_dt = config.start_time
    t_end_dt = config.end_time

    # Create 4 subplots arranged in a vertical stack (4 rows, 1 column)
    fig, axes = plt.subplots(4, 1, figsize=(12, 12))

    # ---------------------------
    # Volumetric strain subplot
    # ---------------------------
    axes[0].plot(t_tide, StressVol, color='black', linestyle='-', linewidth=0.5, label='Volumetric strain')
    axes[0].plot(t_decluster, volumetric_strain, 'ro', markersize=2)
    axes[0].set_ylabel('Volumetric strain', fontsize=16)
    axes[0].grid()  # Enable grid lines for better readability

    # ---------------------------
    # Normal stress subplot
    # ---------------------------
    axes[1].plot(t_tide, StressN, color='black', linestyle='-', linewidth=0.5, label='Normal stress')
    axes[1].plot(t_decluster, normal_stress, 'ro', markersize=2)
    axes[1].set_ylabel('Normal stress (Pa)', fontsize=16)
    axes[1].grid()

    # ---------------------------
    # Shear stress subplot
    # ---------------------------
    axes[2].plot(t_tide, StressS, color='black', linestyle='-', linewidth=0.5, label='Shear stress')
    axes[2].plot(t_decluster, shear_stress, 'ro', markersize=2)
    axes[2].set_ylabel('Shear stress (Pa)', fontsize=16)
    axes[2].grid()

    # ---------------------------
    # Coulomb stress subplot
    # ---------------------------
    axes[3].plot(t_tide, StressCFS, color='black', linestyle='-', linewidth=0.5, label='Coulomb stress')
    axes[3].plot(t_decluster, stress_cfs, 'ro', markersize=2)
    axes[3].set_ylabel('Coulomb stress (Pa)', fontsize=16)
    axes[3].grid()

    # ---------------------------
    # Format the x-axis for all subplots
    # ---------------------------
    # Only display x-axis labels and ticks on the bottom subplot.
    for i, ax in enumerate(axes):
        ax.set_xlim([t_start_dt, t_end_dt])  # Set x-axis limits based on configuration times
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))  # Format x-axis ticks as years
        ax.tick_params(axis='y', labelsize=16)
        if i < len(axes) - 1:
            # Hide x-axis labels for all subplots except the last one
            ax.tick_params(axis='x', which='both', labelbottom=False)
        else:
            # For the bottom subplot, rotate x-axis labels for better readability and add an x-axis label
            ax.tick_params(axis='x', rotation=45, labelsize=16)
            ax.set_xlabel('Time', fontsize=16)

    plt.tight_layout()  # Adjust subplots to fit in the figure area without overlapping
    plt.savefig(config.TM_all_region_stress, format='pdf', dpi=300)  # Save the figure as a PDF file
    # plt.show()  # Uncomment this line to display the plot interactively



def plot_filled_grids_a_value(grid_centers, bin_size, stress_name, data_decluster_rep,ax):
        """
        Plots filled grid cells with a color mapping based on a modulation parameter.

        Parameters:
            grid_centers (ndarray): An (N,5) array with each row:
                [lon_center, lat_center, a_estimated, schuster_value, number_eq].
            bin_size (float): The size of each grid cell (in degrees).
            stress_name (str): Stress name used for the plot title.
            data_decluster_rep (ndarray): Earthquake event data used to plot background events.
                                          It is assumed that column 8 (index 7) contains longitudes
                                          and column 7 (index 6) contains latitudes.
        """
        # Create a new figure and axis
        # fig, ax = plt.subplots()

        # Extract the modulation parameter values from column index 2.
        a_values = grid_centers[:, 2]

        # The MATLAB code sets the color mapping range explicitly.
        pm_min = min(a_values)
        pm_max =  max(a_values)

        # Normalize the modulation values to the interval [0,1]
        norm_values = (a_values - pm_min) / (pm_max - pm_min)

        # Define a colormap with 64 discrete levels.
        cmap = cm.get_cmap('jet', 64)

        # Loop through each grid cell
        for idx in range(grid_centers.shape[0]):
            lon_center = grid_centers[idx, 0]
            lat_center = grid_centers[idx, 1]
            pm_value = grid_centers[idx, 2]
            number_eq = grid_centers[idx, 4]  # For an (N,5) array, the number of events is in column index 4.

            # Only plot cells with at least one event.
            if number_eq == 0:
                continue

            # Define the four vertices of the grid cell
            half_bin = bin_size / 2
            lon_vertices = [lon_center - half_bin, lon_center + half_bin,
                            lon_center + half_bin, lon_center - half_bin]
            lat_vertices = [lat_center - half_bin, lat_center - half_bin,
                            lat_center + half_bin, lat_center + half_bin]

            # Determine the color index based on the normalized value.
            # MATLAB does: round(norm * 63) + 1 (indices 1 to 64). Here we map to an integer 0–63.
            color_idx = int(np.round(norm_values[idx] * 63))
            color_idx = np.clip(color_idx, 0, 63)
            # Use the discrete colormap; note that cmap accepts a float between 0 and 1,
            # so we convert the integer index back to [0,1] by dividing by 63.
            color = cmap(color_idx / 63)

            # Fill the grid cell with the chosen color and no edge.
            ax.fill(lon_vertices, lat_vertices, color=color, edgecolor='none')

        # Set up a ScalarMappable so that the colorbar reflects the same normalization.
        norm = colors.Normalize(vmin=pm_min, vmax=pm_max)
        mappable = cm.ScalarMappable(norm=norm, cmap='jet')
        mappable.set_array([])  # required for older versions of Matplotlib
        cbar = plt.colorbar(mappable, ax=ax)
        cbar.set_label(r'$\alpha$ value', fontsize=18)

        # Plot background earthquake events.
        # MATLAB uses: plot(data_decluster_rep(:,8), data_decluster_rep(:,7), '.k');
        # which corresponds to Python indices 7 (longitude) and 6 (latitude).
        ax.plot(data_decluster_rep[:, 7], data_decluster_rep[:, 6], '.k')

        # Set axis labels and plot title.
        ax.set_xlabel('Longitude', fontsize=18)
        ax.set_ylabel('Latitude', fontsize=18)
        ax.set_title(f'Tidal sensitivity ({stress_name})', fontsize=18)

        # Enhance the plot: grid, tick label size, etc.
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=18)

        # plt.show()



def plot_filled_grids_p_value(grid_centers, bin_size, stress_name, data_decluster_rep,ax):
    """
    Plots filled grid cells based on Schuster test values.

    Parameters:
        grid_centers (ndarray): An (N, 5) array where each row is
            [lon_center, lat_center, a_estimated, schuster_value, number_eq].
        bin_size (float): The grid cell size.
        stress_name (str): Stress name used for setting the plot title.
        data_decluster_rep (ndarray): Earthquake event data for background plotting.
            Assumes that column 8 (index 7) holds event longitudes and column 7 (index 6) holds latitudes.
    """
    # Create a new figure and axis
    # fig, ax = plt.subplots()
    # ax.set_aspect('equal')
    # plt.tight_layout()

    # Extract Schuster values (from column index 3)
    schuster_values = grid_centers[:, 3]

    # Define a two-color colormap:
    # First color (for schuster >= 0.05): gray ([0.8, 0.8, 0.8])
    # Second color (for schuster < 0.05): red ([1, 0, 0])
    cmap = ListedColormap([[0.8, 0.8, 0.8], [1, 0, 0]])

    # Create an index array: if schuster value < 0.05, assign 1 (red), else 0 (gray)
    color_indices = (schuster_values < 0.05).astype(int)

    # Loop through each grid cell
    for idx in range(grid_centers.shape[0]):
        lon_center = grid_centers[idx, 0]
        lat_center = grid_centers[idx, 1]
        # Number of events is assumed to be in the 5th column (index 4)
        number_eq = grid_centers[idx, 4]

        # Skip cells with no events
        if number_eq == 0:
            continue

        # Define the vertices of the grid cell (a square centered at the grid center)
        half_bin = bin_size / 2.0
        lon_vertices = [lon_center - half_bin, lon_center + half_bin,
                        lon_center + half_bin, lon_center - half_bin]
        lat_vertices = [lat_center - half_bin, lat_center - half_bin,
                        lat_center + half_bin, lat_center + half_bin]

        # Choose the color based on the condition:
        # In MATLAB: fill(..., cmap(color_idx(index)+1, :), ...).
        # Here, if color_indices[idx] == 1 (i.e. schuster < 0.05) we pick red; otherwise gray.
        cell_color = cmap.colors[color_indices[idx]]
        ax.fill(lon_vertices, lat_vertices, color=cell_color, edgecolor='none')

    # Plot background earthquake events.
    # MATLAB uses: plot(data_decluster_rep(:,8), data_decluster_rep(:,7), '.k');
    # so in Python, we use indices 7 (longitude) and 6 (latitude).
    ax.plot(data_decluster_rep[:, 7], data_decluster_rep[:, 6], '.k')

    # Create a ScalarMappable to generate a colorbar using our two-color cmap.
    sm = ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, ticks=[0, 1])

    # Set the tick labels exactly as in the MATLAB code:
    # MATLAB: {'Schuster < 0.05', 'Schuster ≥ 0.05'}
    cbar.ax.set_yticklabels(['Schuster < 0.05', 'Schuster ≥ 0.05'])

    # Hide the colorbar (default hidden in MATLAB)
    cbar.ax.set_visible(False)

    # Set axis labels and title
    ax.set_xlabel('Longitude', fontsize=18)
    ax.set_ylabel('Latitude', fontsize=18)
    ax.set_title(f'Schuster Test Value ({stress_name})', fontsize=18)

    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=18)

    # plt.show()




def plot_filled_grids_pha_value(grid_centers, bin_size, stress_name, data_decluster_rep,ax):
    """
    Plots filled grid cells with colors representing the phase shift value.

    Parameters:
        grid_centers (ndarray): An (N, 7) array where each row is
            [lon_center, lat_center, a_estimated, schuster_value, <unused>, pm_value, number_eq].
        bin_size (float): The size of each grid cell.
        stress_name (str): Stress name for setting the plot title.
        data_decluster_rep (ndarray): Background earthquake event data.
            Assumes that column 8 (index 7) holds event longitudes and column 7 (index 6) holds latitudes.
    """
    # Create a new figure and get the current axis.
    # fig, ax = plt.subplots()

    # Extract the phase shift values (pm_value) from column index 5.
    pm_values = grid_centers[:, 5]
    pm_min = np.min(pm_values)
    pm_max = np.max(pm_values)

    # Loop through each grid cell.
    for i in range(grid_centers.shape[0]):
        lon_center = grid_centers[i, 0]
        lat_center = grid_centers[i, 1]
        pm_value = grid_centers[i, 5]  # phase shift value
        number_eq = grid_centers[i, 6]  # number of earthquakes in this grid cell

        # Only plot cells with at least one earthquake.
        if number_eq == 0:
            continue

        # Determine the vertices of the grid cell (a square centered at (lon_center, lat_center)).
        half_bin = bin_size / 2.0
        lon_vertices = [lon_center - half_bin, lon_center + half_bin,
                        lon_center + half_bin, lon_center - half_bin]
        lat_vertices = [lat_center - half_bin, lat_center - half_bin,
                        lat_center + half_bin, lat_center + half_bin]

        # Normalize the phase shift value.
        if pm_max != pm_min:
            norm_value = (pm_value - pm_min) / (pm_max - pm_min)
        else:
            norm_value = 0

        # Compute a color that transitions from red (low pm_value) to blue (high pm_value).
        # MATLAB code uses: color = [1 - norm_value, 0, norm_value]
        color = [1 - norm_value, 0, norm_value]

        # Fill the grid cell with the computed color (no edge color).
        ax.fill(lon_vertices, lat_vertices, color=color, edgecolor='none')

    # Set up a colorbar using a jet colormap with the same data range.
    norm = mcolors.Normalize(vmin=pm_min, vmax=pm_max)
    sm = ScalarMappable(cmap=plt.cm.jet, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Phase Angle (degrees)', fontsize=18)

    # Plot the background earthquake events.
    # MATLAB: plot(data_decluster_rep(:,8), data_decluster_rep(:,7), '.k');
    ax.plot(data_decluster_rep[:, 7], data_decluster_rep[:, 6], '.k')

    # Set axis labels and plot title.
    ax.set_xlabel('Longitude', fontsize=18)
    ax.set_ylabel('Latitude', fontsize=18)
    ax.set_title(f'Phase shift ({stress_name})', fontsize=18)

    # Enhance the plot appearance.
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=18)

    # plt.show()


def plot_filled_grids_amp_value(grid_centers, bin_size, stress_name, data_decluster_rep,ax):
    """
    Plots filled grid cells with colors representing the amplitude (a_estimated) from grid_centers.

    Parameters:
        grid_centers (ndarray): An (N,5) array where each row is
            [lon_center, lat_center, a_estimated, schuster_value, number_eq].
        bin_size (float): The grid cell size.
        stress_name (str): Stress name for setting the plot title.
        data_decluster_rep (ndarray): Data for background earthquake events.
            Assumes that column 8 (Python index 7) holds event longitudes and column 7 (Python index 6) holds latitudes.
    """

    # Extract amplitude values (a_estimated) from column index 2
    amp_values = grid_centers[:, 4]
    amp_min = np.min(amp_values)
    amp_max = np.max(amp_values)

    # Normalize amplitude values to [0, 1]
    norm_values = (amp_values - amp_min) / (amp_max - amp_min) if amp_max != amp_min else np.zeros_like(amp_values)

    # Create a discrete jet colormap with 64 levels
    cmap = plt.get_cmap('jet', 64)

    # Loop over each grid cell
    for i in range(grid_centers.shape[0]):
        lon_center = grid_centers[i, 0]
        lat_center = grid_centers[i, 1]
        amp_value = grid_centers[i, 2]
        # Number of events is in column index 4
        number_eq = grid_centers[i, 4]

        # Only plot cells with at least one earthquake
        if number_eq == 0:
            continue

        # Determine the four vertices of the grid cell (square)
        half_bin = bin_size / 2.0
        lon_vertices = [lon_center - half_bin, lon_center + half_bin,
                        lon_center + half_bin, lon_center - half_bin]
        lat_vertices = [lat_center - half_bin, lat_center - half_bin,
                        lat_center + half_bin, lat_center + half_bin]

        # Normalize the amplitude value and convert to a discrete color index
        norm_val = norm_values[i]
        color_idx = int(round(norm_val * 63))
        color_idx = np.clip(color_idx, 0, 63)
        # Obtain the color from the colormap
        color = cmap(color_idx / 63.0)

        # Fill the grid cell with the chosen color (no edge)
        ax.fill(lon_vertices, lat_vertices, color=color, edgecolor='none')

    # Set up a ScalarMappable for the colorbar using the same jet colormap and data range
    norm = mcolors.Normalize(vmin=amp_min, vmax=amp_max)
    sm = cm.ScalarMappable(norm=norm, cmap='jet')
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('pm value', fontsize=18)

    # Plot background earthquake events
    # MATLAB: plot(data_decluster_rep(:,8), data_decluster_rep(:,7), '.k')
    ax.plot(data_decluster_rep[:, 7], data_decluster_rep[:, 6], '.k')

    # Set axis labels and title
    ax.set_xlabel('Longitude', fontsize=18)
    ax.set_ylabel('Latitude', fontsize=18)
    ax.set_title(f'Amplitude ({stress_name})', fontsize=18)

    # Enhance the plot appearance
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.grid(True)
    plt.box(on=True)

