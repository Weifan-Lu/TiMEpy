import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime

import b_positive as bp
import matplotlib as mpl
import visualization as vis

# Set the global font size to 12
mpl.rcParams['font.size'] = 12


def ana4_b_pos(config):
    plt.close('all')

    # 1. Load data
    data = np.loadtxt(config.data_select)

    # 2. Sort by time (assuming the first six columns are Year, Month, Day, Hour, Minute, Second)
    sort_indices = np.lexsort((data[:, 5], data[:, 4], data[:, 3],
                               data[:, 2], data[:, 1], data[:, 0]))
    sorted_data = data[sort_indices]

    # 3. Extract magnitude data (assuming magnitude is in the 10th column, i.e., index 9)
    magnitudes = sorted_data[:, 9]

    # 4. Convert time data to datetime objects (assuming the first six columns are Year, Month, Day, Hour, Minute, Second)
    t_data = np.array([
        datetime.datetime(int(row[0]), int(row[1]), int(row[2]),
                          int(row[3]), int(row[4]), int(row[5]))
        for row in sorted_data
    ])

    # 5. Get the time range and moving window parameters
    dt_start = config.start_time  # Should be a datetime object
    dt_end = config.end_time  # Should be a datetime object
    tmin = mdates.date2num(dt_start)
    tmax = mdates.date2num(dt_end)
    step_days = config.step_days  # Step size (unit: days)
    window_days = config.window_days  # Window length (unit: days)

    # Lists to store results
    window_centers = []
    b_values = []
    b_std_all = []

    # 6. Sliding window traversal
    # Generate t_start as a matplotlib date number using np.arange
    t_data_num = mdates.date2num(t_data)  # For easy comparison
    for t_start in np.arange(tmin, tmax - window_days + 1e-6, step_days):
        t_end = t_start + window_days
        # Select observed data within the current window
        indices = np.where((t_data_num >= t_start) & (t_data_num <= t_end))[0]
        mag_cut = magnitudes[indices]

        # If there is enough data in the window, compute the b-value; otherwise, set it to NaN
        if len(mag_cut) > 0:
            b_val = bp.b_positive_discrete(mag_cut, delta=0.05)
            b_mean, b_std = bp.b_positive_discrete_bootstrap(magnitudes, delta=0.05, n_boot=1000)
        else:
            b_val = np.nan

        # Save the center time of the window (using the midpoint of the start and end times)
        window_center = t_start + window_days
        window_centers.append(window_center)
        b_values.append(b_val)
        b_std_all.append(b_std)

    # Convert window centers from numerical format to datetime format
    window_dates = mdates.num2date(window_centers)
    komo_t = config.mainshock_t

    # Create two subplots sharing the x-axis
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=None, figsize=(8, 6))

    # Subplot 1: Plot b-value variations over time with error bars
    ax1.errorbar(
        window_dates, b_values, yerr=b_std_all, fmt='o',
        markersize=6, markerfacecolor='royalblue', markeredgecolor='black',  # Enhance marker appearance
        capsize=4, capthick=1.2, elinewidth=2, alpha=0.8, color='royalblue',  # Enhance error bars
        linestyle='-')

    # Set x-axis to display dates
    # ax1.xaxis_date()
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    # fig.autofmt_xdate()  # Auto-rotate date labels

    ax1.set_ylabel('b-positive')
    ax1.axvline(komo_t, color='r', linewidth=2)
    vis.set_common_elements(ax1, window_centers, 'Time', 'b-positive', 'Time vs. b-positive', font_size=14)

    plt.tight_layout()
    plt.savefig(config.TM_earthquake_b_positive, format='pdf')
    plt.show()
