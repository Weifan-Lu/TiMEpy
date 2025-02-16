import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from scipy.ndimage import uniform_filter

def filtered_catalog_reg(data, lat_center, lon_center, side_length_x, side_length_y, angle_deg,
                         depth_cut=None, mag_cut=None, start_time=None, end_time=None):
    """
    Filters earthquake catalog data within a rotated rectangle.

    Parameters:
        data         : NumPy array with columns
                       [year, month, day, hour, minute, second, latitude, longitude, depth, magnitude, ...]
        lat_center   : Center latitude of the study region.
        lon_center   : Center longitude of the study region.
        side_length_x: Length of the rectangle in longitude (degrees).
        side_length_y: Length of the rectangle in latitude (degrees).
        angle_deg    : Rotation angle in degrees.
        depth_cut    : Maximum depth (if provided).
        mag_cut      : Minimum magnitude (if provided).
        start_time   : Earliest event time (if provided; either a datetime or ISO string).
        end_time     : Latest event time (if provided; either a datetime or ISO string).

    Returns:
        filtered_data: Filtered data array.
        rect_vertices: A 4x2 array with the rectangle's vertices as [latitude, longitude].
    """
    # Compute half sizes
    half_width = side_length_x / 2.0
    half_height = side_length_y / 2.0

    # Define rectangle corners before rotation (relative offsets: [dx, dy])
    corners = np.array([
        [-half_width, -half_height],
        [half_width, -half_height],
        [half_width, half_height],
        [-half_width, half_height]
    ])

    # Rotation matrix (angle in radians)
    angle_rad = np.deg2rad(angle_deg)
    R = np.array([
        [np.cos(angle_rad), -np.sin(angle_rad)],
        [np.sin(angle_rad), np.cos(angle_rad)]
    ])
    # Rotate corners
    rotated_corners = (R @ corners.T).T
    # Calculate rectangle vertices in geographic coordinates:
    # latitude = lat_center + dy, longitude = lon_center + dx
    rect_vertices = np.column_stack((lat_center + rotated_corners[:, 1],
                                     lon_center + rotated_corners[:, 0]))

    # Preliminary filtering using the envelope of the rotated rectangle
    lat_env_min = lat_center + np.min(rotated_corners[:, 1])
    lat_env_max = lat_center + np.max(rotated_corners[:, 1])
    lon_env_min = lon_center + np.min(rotated_corners[:, 0])
    lon_env_max = lon_center + np.max(rotated_corners[:, 0])

    # Get columns (Python uses 0-indexing)
    lat = data[:, 6]
    lon = data[:, 7]
    depth = data[:, 8]
    mag = data[:, 9]

    # Convert date/time columns (columns 0-5) to datetime objects
    times = np.array([datetime(int(row[0]), int(row[1]), int(row[2]),
                               int(row[3]), int(row[4]), int(row[5]))
                      for row in data])

    # Apply spatial envelope filter
    idx = ((lat >= lat_env_min) & (lat <= lat_env_max) &
           (lon >= lon_env_min) & (lon <= lon_env_max))

    # If additional filters are provided, apply them here
    if depth_cut is not None:
        idx &= (depth <= depth_cut)
    if mag_cut is not None:
        idx &= (mag >= mag_cut)
    if (start_time is not None) and (end_time is not None):
        # Ensure start_time and end_time are datetime objects
        if not isinstance(start_time, datetime):
            start_time = datetime.fromisoformat(start_time)
        if not isinstance(end_time, datetime):
            end_time = datetime.fromisoformat(end_time)
        idx &= ((times >= start_time) & (times <= end_time))

    data_pre = data[idx, :]
    times_pre = times[idx]

    # Now do the precise filtering: rotate the coordinates back
    rel_lon = data_pre[:, 7] - lon_center
    rel_lat = data_pre[:, 6] - lat_center
    rel_coords = np.column_stack((rel_lon, rel_lat))

    # Inverse rotation (transpose of R)
    inv_R = R.T
    rotated_coords = rel_coords @ inv_R  # shape (n,2)

    # Select events within the unrotated rectangle boundaries
    idx2 = ((rotated_coords[:, 0] >= -half_width) & (rotated_coords[:, 0] <= half_width) &
            (rotated_coords[:, 1] >= -half_height) & (rotated_coords[:, 1] <= half_height))

    filtered_data = data_pre[idx2, :]

    return filtered_data, rect_vertices


