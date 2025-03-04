import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as patches
from datetime import datetime
import main_select_catalog as msc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from obspy.imaging.beachball import beach

# Helper Function: Calculate the Position of the Count Text in the Plot
def calc_text_position(lon_vals, lat_vals):
    text_lon = np.min(lon_vals) - 0.03 * (np.max(lon_vals) - np.min(lon_vals))
    text_lat = np.max(lat_vals) + 0.08 * (np.max(lat_vals) - np.min(lat_vals))
    return text_lon, text_lat


def pre1_select_catalog(config):
    """
    Functionality:
      1. Load and sort the raw data by time.
      2. Perform stepwise filtering of the data based on spatial constraints, depth, magnitude, and time.
      3. Output the filtered catalog file and generate plots of the raw and filtered data.

    Parameters (the config object must contain the following attributes):
      - data_original: Path to the raw data file (text format, numerical data).
      - lat_center, lon_center: Latitude and longitude of the center point.
      - side_length_x, side_length_y: Side lengths of the rectangle in the x and y directions.
      - angle_deg: Rotation angle of the rectangle (in degrees).
      - start_time, end_time: Time range (ISO format string or datetime object).
      - depth_cut: Depth filtering threshold.
      - data_select: Output filename for the filtered data.
      - select_data_fig: Filename for saving the plot.
    """
    # Close all existing plots

    plt.close('all')

    # ----- 1. Load and sort raw data -----
    # Assume the data file contains the following columns per row:
    # year, month, day, hour, minute, sec, lat, lon, depth, mag, ...

    data_original = np.loadtxt(config.data_original)

    # Sort based on the first six columns (time)
    sort_keys = (data_original[:, 5], data_original[:, 4], data_original[:, 3],
                 data_original[:, 2], data_original[:, 1], data_original[:, 0])
    sorted_idx = np.lexsort(sort_keys)
    data_original = data_original[sorted_idx]

    # ----- 2. Get Parameters-----
    lat_center = config.lat_center
    lon_center = config.lon_center
    side_length_x = config.side_length_x
    side_length_y = config.side_length_y
    angle_deg = config.angle_deg
    start_time = config.start_time
    end_time = config.end_time
    depth_cut = config.depth_cut

    # ----- 3. First Step Filtering: Perform Spatial Filtering on the Entire Rectangular Area-----
    filtered_data_original, _ = msc.filtered_catalog_reg(
        data_original, lat_center, lon_center, side_length_x*2, side_length_y*2, angle_deg
    )

    # ----- 5. Third Step Filtering: Further Filter by Depth, Magnitude, and Time-----
    mag_cut = -10  # Reset the Magnitude Filtering Parameters
    filtered_data_mc, rect_vertices = msc.filtered_catalog_reg(
        filtered_data_original, lat_center, lon_center,
        side_length_x, side_length_y, angle_deg,
        depth_cut, mag_cut, start_time, end_time
    )

    # Convert the Time of Each Filtered Record (First 6 Columns) to a datetime Object
    filtered_event_time = np.array([
        datetime(int(row[0]), int(row[1]), int(row[2]), int(row[3]), int(row[4]), int(row[5]))
        for row in filtered_data_mc
    ])

    num_original = filtered_data_original.shape[0]
    num_filtered = filtered_data_mc.shape[0]

    # ----- 6. Write the Filtered Data to a File-----
    with open(config.data_select, 'w') as fid:
        for row in filtered_data_mc:
            # output：year month day hour minute sec lat lon depth mag
            line = (
                f"{int(round(row[0]))} {int(round(row[1])):02d} {int(round(row[2])):02d} "
                f"{int(round(row[3])):02d} {int(round(row[4])):02d} {row[5]:05.2f} "
                f"{row[6]} {row[7]} {row[8]:04.2f} {row[9]:04.2f}\n"
            )
            fid.write(line)

    # ----- 7. Plotting -----
    # Define font size parameters
    FONT_LABEL = 16  # Font size for axis labels
    FONT_TITLE = 16  # Font size for the chart title
    FONT_TEXT = 16  # Font size for text displayed in the plot
    COLORBAR_LABEL = 12  # Font size for the colorbar label

    # Calculate the rectangular envelope (used to display the filtered region in the plot)

    #  Create a 2x2 subplot
    fig, axs = plt.subplots(2, 2, figsize=(12, 12))

    # --- Subplot 1: Original data spatial distribution (colored by depth) ---
    ax = axs[0, 0]
    sc1 = ax.scatter(
        filtered_data_original[:, 7], filtered_data_original[:, 6],
        s=17, c=filtered_data_original[:, 8], cmap='jet', alpha=0.7
    )
    ax.set_xlabel('Longitude', fontsize=FONT_LABEL)
    ax.set_ylabel('Latitude', fontsize=FONT_LABEL)
    ax.set_title('Background data', fontsize=FONT_TITLE)
    ax.axis('equal')
    ax.grid(True)
    cbar1 = fig.colorbar(sc1, ax=ax)
    cbar1.set_label('Depth (km)', fontsize=COLORBAR_LABEL)
    text_lon, text_lat = calc_text_position(filtered_data_original[:, 7],
                                            filtered_data_original[:, 6])
    ax.text(text_lon, text_lat, f'Count: {num_original}', fontsize=FONT_TEXT,
            fontweight='bold', color='k',
            bbox=dict(facecolor='w', edgecolor='k'))

    # --- Subplot 2: Filtered data (original data in gray, selected events colored by depth) ---
    ax = axs[0, 1]
    ax.scatter(filtered_data_original[:, 7], filtered_data_original[:, 6],
               s=10, color='gray', alpha=0.7)
    sc2 = ax.scatter(
        filtered_data_mc[:, 7], filtered_data_mc[:, 6],
        s=10, c=filtered_data_mc[:, 8], cmap='jet', alpha=0.7
    )
    ax.set_xlabel('Longitude', fontsize=FONT_LABEL)
    ax.set_ylabel('Latitude', fontsize=FONT_LABEL)
    ax.set_title('Selected data', fontsize=FONT_TITLE)
    ax.axis('equal')
    ax.grid(True)
    cbar2 = fig.colorbar(sc2, ax=ax)
    cbar2.set_label('Depth (km)', fontsize=COLORBAR_LABEL)
    text_lon, text_lat = calc_text_position(filtered_data_mc[:, 7],
                                            filtered_data_mc[:, 6])
    ax.text(text_lon, text_lat, f'Count: {num_filtered}', fontsize=FONT_TEXT,
            fontweight='bold', color='k',
            bbox=dict(facecolor='w', edgecolor='k'))

    # --- Subplot 3: Time-Magnitude Scatter Plot ---

    ax = axs[1, 0]
    ax.scatter(filtered_event_time, filtered_data_mc[:, 9],
               s=30, color='r', alpha=0.7)
    ax.set_xlabel('Time', fontsize=FONT_LABEL)
    ax.set_ylabel('Magnitude', fontsize=FONT_LABEL)
    ax.set_title('Magnitude vs. Time', fontsize=FONT_TITLE)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.grid(True)

    # --- Subplot 4: Cumulative Number of Earthquakes Over Time ---

    ax = axs[1, 1]
    cum_count = np.arange(1, len(filtered_event_time) + 1)
    ax.plot(filtered_event_time, cum_count, '-r', linewidth=2)
    ax.set_xlabel('Time', fontsize=FONT_LABEL)
    ax.set_ylabel('Cumulative Count', fontsize=FONT_LABEL)
    ax.set_title('Cumulative Earthquake Count', fontsize=FONT_TITLE)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.grid(True)

    plt.tight_layout()
    fig.savefig(config.select_data_fig, dpi=150, bbox_inches='tight')
    plt.show()

    #
    #
    # # Create the figure: Map subplots on the left, earthquake activity plots on the right
    # fig = plt.figure(figsize=(18, 6))
    # ax_map = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
    # ax_mt = fig.add_subplot(1, 2, 2)
    #
    # # --------------------- Map Subplots ---------------------
    # # Add map background elements
    #
    # ax_map.add_feature(cfeature.LAND, facecolor='lightgray')
    # ax_map.add_feature(cfeature.OCEAN, facecolor='lightblue')
    # ax_map.add_feature(cfeature.COASTLINE, linewidth=1)
    # ax_map.add_feature(cfeature.BORDERS, linestyle=':')
    #
    # # Plot scatter plot of original data (specify transform to ensure correct projection)
    # sc = ax_map.scatter(
    #     data_original[:, 7], data_original[:, 6],
    #     s=5, c='grey',
    #     transform=ccrs.PlateCarree(), rasterized=True
    # )
    # # Plot scatter plot of original data (specify transform to ensure correct projection)
    # sc = ax_map.scatter(
    #     filtered_data_mc[:, 7], filtered_data_mc[:, 6],
    #     s=5, c=filtered_data_mc[:, 8], cmap='jet', alpha=0.7,
    #     transform=ccrs.PlateCarree(), rasterized=True
    # )
    #
    # # Draw a rectangle representing the filtered area
    # rect = patches.Rectangle(
    #     (lon_min, lat_min), width, height,
    #     linewidth=2, edgecolor='r', facecolor='none',
    #     transform=ccrs.PlateCarree()
    # )
    # ax_map.add_patch(rect)
    #
    # # Set coordinate labels
    # ax_map.set_xlabel('Longitude', fontsize=FONT_LABEL)
    # ax_map.set_ylabel('Latitude', fontsize=FONT_LABEL)
    #
    # # Set map display range
    # bin_map_s = 0.5
    # ax_map.set_extent(
    #     [lon_min - bin_map_s, lon_min + width + bin_map_s, lat_min - bin_map_s, lat_min + height + bin_map_s],
    #     crs=ccrs.PlateCarree()
    # )
    #
    # cax = inset_axes(ax_map,
    #                  width="2%",  # Colorbar width, 3% of ax_map width
    #                  height="30%",  # Colorbar height, 50% of ax_map height
    #                  loc='lower left',
    #                  bbox_to_anchor=(0.85, 0.05, 1, 1),  # Adjust colorbar position in ax_map
    #                  bbox_transform=ax_map.transAxes,
    #                  borderpad=0)
    #
    # # Draw the colorbar on cax
    # cbar = fig.colorbar(sc, cax=cax, orientation='vertical')
    # cbar.set_label('Depth (km)', fontsize=COLORBAR_LABEL)
    # # Set the background of the colorbar to white
    # cbar.ax.set_facecolor('white')
    # # Add a small inset map in the top right corner
    # ax_inset = fig.add_axes([0.28, 0.68, 0.2, 0.3], projection=ccrs.PlateCarree())
    # # ax_inset.set_global()  # Show global map
    # ax_inset.add_feature(cfeature.OCEAN, facecolor='lightblue', zorder=0, rasterized=True)
    # ax_inset.add_feature(cfeature.LAND, facecolor='lightgray', zorder=1, rasterized=True)
    # ax_inset.add_feature(cfeature.COASTLINE, linewidth=1, zorder=2, rasterized=True)
    # ax_inset.add_feature(cfeature.BORDERS, linestyle=':', zorder=2, rasterized=True)
    # ax_inset.add_feature(cfeature.STATES, edgecolor='black', zorder=3, rasterized=True)
    #
    # # Draw a red rectangle on the inset map to mark the study area
    # rect_inset = patches.Rectangle(
    #     (lon_min, lat_min), width, height,
    #     linewidth=1.5, edgecolor='r', facecolor='none',
    #     transform=ccrs.PlateCarree()
    # )
    # ax_inset.add_patch(rect_inset)
    # # Adjust inset map display range to make the study area more visible
    # ax_inset.set_extent(
    #     [lon_min - 5, lon_min + width + 5, lat_min - 5, lat_min + height + 5],
    #     crs=ccrs.PlateCarree()
    # )
    #
    # strike = config.str
    # dip = config.dip
    # rake = config.slip
    # main_lat = config.main_lat
    # main_lon = config.main_lon
    # width_ball = 0.1
    # bball = beach([strike, dip, rake],
    #               xy=(main_lon, main_lat),
    #               width=width_ball,
    #               linewidth=1,
    #               facecolor='red')
    # # Key: Set the coordinates of the beachball to geographic coordinates (PlateCarree)
    # bball.set_transform(ccrs.PlateCarree())
    #
    # # 6. Add beachball to ax_map
    # ax_map.add_collection(bball)
    #
    # gl = ax_map.gridlines(draw_labels=True, linestyle="--", linewidth=0.5, color="gray")
    # gl.xlocator = mticker.MultipleLocator(0.4)
    # gl.ylocator = mticker.MultipleLocator(0.4)
    # gl.right_labels = False
    # gl.top_labels = False
    #
    # # --------------------- Seismic Activity Plot (Right Panel) ---------------------
    # # ax_mt.set_title("Magnitude & Cumulative Count vs Time", fontsize=14)
    # ax_mt.vlines(filtered_event_time, ymin=-2, ymax=filtered_data_mc[:, 9], colors='gray', linestyles='-',
    #              linewidth=0.5)
    # ax_mt.plot(filtered_event_time, filtered_data_mc[:, 9], 'o', color='black', markersize=1, label='Magnitude')
    #
    # ax_mt.set_xlabel("Time")
    # ax_mt.set_ylabel("Magnitude", color='black')
    # ax_mt.tick_params(axis='y', labelcolor='black')
    #
    # # Format time axis
    # ax_mt.xaxis.set_major_locator(mdates.AutoDateLocator())
    # ax_mt.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax_mt.xaxis.get_major_locator()))
    # fig.autofmt_xdate()
    #
    # # Plot cumulative earthquake count
    # cumulative_counts = np.arange(1, len(filtered_event_time) + 1)
    # ax_mt2 = ax_mt.twinx()
    # ax_mt2.plot(filtered_event_time, cumulative_counts, 'o-', color='blue', markersize=2, label='Cumulative')
    # ax_mt2.set_ylabel("Cumulative Count", color='blue')
    # ax_mt2.tick_params(axis='y', labelcolor='blue')
    #
    # # Merge legends
    # lines1, labels1 = ax_mt.get_legend_handles_labels()
    # lines2, labels2 = ax_mt2.get_legend_handles_labels()
    # ax_mt2.legend(lines1 + lines2, labels1 + labels2, loc='best')
    #
    # # Add subplot annotations (a) and (b) to their respective top-left corners
    # ax_map.text(-0.1, 0.98, '(a)', transform=ax_map.transAxes,
    #             fontsize=18, verticalalignment='top')
    # ax_mt.text(-0.07, 0.98, '(b)', transform=ax_mt.transAxes,
    #            fontsize=18, verticalalignment='top')
    # plt.tight_layout()
    # fig.savefig(config.select_data_fig_map, format='pdf', dpi=300)
    # plt.show()
