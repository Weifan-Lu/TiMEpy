import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as patches
from datetime import datetime
import main_select_catalog as msc


# 辅助函数：计算图中显示计数文本的位置
def calc_text_position(lon_vals, lat_vals):
    text_lon = np.min(lon_vals) - 0.03 * (np.max(lon_vals) - np.min(lon_vals))
    text_lat = np.max(lat_vals) + 0.08 * (np.max(lat_vals) - np.min(lat_vals))
    return text_lon, text_lat


def pre1_select_catalog_fault(config):
    """
    MATLAB函数 pre1_select_catalog_fault 的 Python 翻译版。

    功能：
      1. 加载并按时间排序原始数据
      2. 分步对数据进行空间、深度、震级及时间过滤
      3. 输出过滤后的目录文件，并绘制原始和过滤后数据的图像

    参数（config 对象需包含以下属性）：
      - data_original: 原始数据文件路径（文本格式，数据为数值型）
      - lat_center, lon_center: 中心点经纬度
      - side_length_x, side_length_y: 矩形的 x、y 方向边长
      - angle_deg: 矩形旋转角度（度）
      - start_time, end_time: 时间范围（ISO格式字符串或 datetime 对象）
      - depth_cut: 深度过滤阈值
      - data_select: 输出的过滤后数据文件名
      - select_data_fig: 保存的图像文件名
    """
    # 关闭所有已有图像
    plt.close('all')

    # ----- 1. 加载并排序原始数据 -----
    # 假设数据文件每行包含：year, month, day, hour, minute, sec, lat, lon, depth, mag, ...
    data_original = np.loadtxt(config.data_original)

    # 根据前六列（时间）排序
    sort_keys = (data_original[:, 5], data_original[:, 4], data_original[:, 3],
                 data_original[:, 2], data_original[:, 1], data_original[:, 0])
    sorted_idx = np.lexsort(sort_keys)
    data_original = data_original[sorted_idx]

    # ----- 2. 获取参数 -----
    lat_center = config.lat_center
    lon_center = config.lon_center
    side_length_x = config.side_length_x
    side_length_y = config.side_length_y
    angle_deg = config.angle_deg
    start_time = config.start_time
    end_time = config.end_time
    depth_cut = config.depth_cut

    # ----- 3. 第一步过滤：对整个矩形区域进行空间过滤 -----
    filtered_data_original, _ = msc.filtered_catalog_reg(
        data_original, lat_center, lon_center, side_length_x, side_length_y, angle_deg
    )

    # ----- 4. 第二步过滤：对半尺寸矩形区域进行过滤 -----
    filtered_data, _ = msc.filtered_catalog_reg(
        filtered_data_original, lat_center, lon_center,
        side_length_x / 2, side_length_y / 2, angle_deg
    )

    # 构建简化后的目录数据：[year, month, lat, lon, depth, mag]
    catamCatalog = np.column_stack((
        filtered_data[:, 0], filtered_data[:, 1],
        filtered_data[:, 6], filtered_data[:, 7],
        filtered_data[:, 8], filtered_data[:, 9]
    ))

    # ----- 5. 第三步过滤：进一步按深度、震级和时间过滤 -----
    mag_cut = -10  # 重新设定震级过滤参数
    filtered_data_mc, rect_vertices = msc.filtered_catalog_reg(
        filtered_data, lat_center, lon_center,
        side_length_x / 2, side_length_y / 2, angle_deg,
        depth_cut, mag_cut, start_time, end_time
    )

    # 将过滤后的每条记录的时间（前6列）转换为 datetime 对象
    filtered_event_time = np.array([
        datetime(int(row[0]), int(row[1]), int(row[2]), int(row[3]), int(row[4]), int(row[5]))
        for row in filtered_data_mc
    ])

    num_original = filtered_data_original.shape[0]
    num_filtered = filtered_data_mc.shape[0]

    # ----- 6. 写入过滤后的数据到文件 -----
    with open(config.data_select, 'w') as fid:
        for row in filtered_data_mc:
            # 格式化输出：year month day hour minute sec lat lon depth mag
            line = (
                f"{int(round(row[0]))} {int(round(row[1])):02d} {int(round(row[2])):02d} "
                f"{int(round(row[3])):02d} {int(round(row[4])):02d} {row[5]:05.2f} "
                f"{row[6]} {row[7]} {row[8]:04.2f} {row[9]:04.2f}\n"
            )
            fid.write(line)

    # ----- 7. 绘图 -----
    # 定义字体大小参数
    FONT_LABEL = 16  # 坐标轴标签字体大小
    FONT_TITLE = 16  # 标题字体大小
    FONT_TEXT = 16  # 图中显示文本的字体大小
    COLORBAR_LABEL = 16  # colorbar标签字体大小

    # 计算矩形包络（用于绘图中显示过滤区域）的参数
    lon_min = np.min(rect_vertices[:, 1])
    lat_min = np.min(rect_vertices[:, 0])
    width = np.max(rect_vertices[:, 1]) - lon_min
    height = np.max(rect_vertices[:, 0]) - lat_min

    # 创建 2x2 子图
    fig, axs = plt.subplots(2, 2, figsize=(12, 12))

    # --- 子图 1：原始数据空间分布（以深度着色） ---
    ax = axs[0, 0]
    sc1 = ax.scatter(
        filtered_data_original[:, 7], filtered_data_original[:, 6],
        s=17, c=filtered_data_original[:, 8], cmap='jet', alpha=0.7
    )
    # 绘制表示过滤区域的矩形
    rect = patches.Rectangle((lon_min, lat_min), width, height,
                             linewidth=2, edgecolor='r', facecolor='none')
    ax.add_patch(rect)
    ax.set_xlabel('Longitude', fontsize=FONT_LABEL)
    ax.set_ylabel('Latitude', fontsize=FONT_LABEL)
    ax.set_title('Original Data', fontsize=FONT_TITLE)
    ax.axis('equal')
    ax.grid(True)
    cbar1 = fig.colorbar(sc1, ax=ax)
    cbar1.set_label('Depth (km)', fontsize=COLORBAR_LABEL)
    text_lon, text_lat = calc_text_position(filtered_data_original[:, 7],
                                            filtered_data_original[:, 6])
    ax.text(text_lon, text_lat, f'Count: {num_original}', fontsize=FONT_TEXT,
            fontweight='bold', color='k',
            bbox=dict(facecolor='w', edgecolor='k'))

    # --- 子图 2：过滤后数据（原始数据为灰色，选中事件以深度着色） ---
    ax = axs[0, 1]
    ax.scatter(filtered_data_original[:, 7], filtered_data_original[:, 6],
               s=10, color='gray', alpha=0.7)
    sc2 = ax.scatter(
        filtered_data_mc[:, 7], filtered_data_mc[:, 6],
        s=10, c=filtered_data_mc[:, 8], cmap='jet', alpha=0.7
    )
    ax.set_xlabel('Longitude', fontsize=FONT_LABEL)
    ax.set_ylabel('Latitude', fontsize=FONT_LABEL)
    ax.set_title('Filtered Data', fontsize=FONT_TITLE)
    ax.axis('equal')
    ax.grid(True)
    cbar2 = fig.colorbar(sc2, ax=ax)
    cbar2.set_label('Depth (km)', fontsize=COLORBAR_LABEL)
    text_lon, text_lat = calc_text_position(filtered_data_mc[:, 7],
                                            filtered_data_mc[:, 6])
    ax.text(text_lon, text_lat, f'Count: {num_filtered}', fontsize=FONT_TEXT,
            fontweight='bold', color='k',
            bbox=dict(facecolor='w', edgecolor='k'))

    # --- 子图 3：时间-震级散点图 ---
    ax = axs[1, 0]
    ax.scatter(filtered_event_time, filtered_data_mc[:, 9],
               s=30, color='r', alpha=0.7)
    ax.set_xlabel('Time', fontsize=FONT_LABEL)
    ax.set_ylabel('Magnitude', fontsize=FONT_LABEL)
    ax.set_title('Magnitude vs. Time', fontsize=FONT_TITLE)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.grid(True)

    # --- 子图 4：累计地震数随时间变化 ---
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
    # plt.show()
