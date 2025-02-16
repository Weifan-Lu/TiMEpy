import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from scipy.ndimage import uniform_filter
import main_declustering_nna as mdnna

def pre2_decluster_NNA(config):
    """
    从 params（例如 input_params）中读取参数与数据文件，
    提取年、月、日、时、分、秒、经纬度、深度、震级，
    并将时间转换为数值型日期（单位：天）。

    参数:
      params: dict 类型，需包含如下键：
            - 'start_time': 开始时间（数值或其他标识，可用于后续计算）
            - 'data_select': 数据文件路径（文本格式，每行包含 10 列数据，
                             顺序为 [year, month, day, hour, minute, sec, lat, lon, dep, mag]）
            - 'df', 'b', 'p', 'q', 'eta0': 算法参数
    返回:
      无返回值（函数内绘图并将去相关目录保存到文件）
    """
    # 读取参数
    start_time = config.start_time
    data_file = config.data_select
    df = config.df
    b = config.b
    p = config.p
    q = config.q
    eta0 = config.eta0

    # 载入数据（假定为文本格式，每行 10 列）
    data = np.loadtxt(data_file)
    results = mdnna.decluster_nna(start_time, data, df, b, p, q, eta0)
    Declustered_Catalog = results['Declustered_Catalog']
    Dataobj = results['Dataobj']
    t = results['t']

    # ======== 将去相关目录输出到文件 ============
    output_filename = config.data_select_decluster
    with open(output_filename, 'w') as fid:
        # 格式为：'%d %02d %02d %02d %02d %05.2f %f %f %04.2f %04.2f %d\n'
        for i in range(Declustered_Catalog.shape[0]):
            yr = int(round(Declustered_Catalog[i, 0]))
            mon = int(round(Declustered_Catalog[i, 1]))
            day_ = int(round(Declustered_Catalog[i, 2]))
            hr = int(round(Declustered_Catalog[i, 3]))
            minu = int(round(Declustered_Catalog[i, 4]))
            sec_ = Declustered_Catalog[i, 5]
            lat_val = Declustered_Catalog[i, 6]
            lon_val = Declustered_Catalog[i, 7]
            dep_val = Declustered_Catalog[i, 8]
            mag_val = Declustered_Catalog[i, 9]
            extra = int(round(Declustered_Catalog[i, 10]))
            line = f"{yr} {mon:02d} {day_:02d} {hr:02d} {minu:02d} {sec_:05.2f} {lat_val:f} {lon_val:f} {dep_val:04.2f} {mag_val:04.2f} {extra:d}\n"
            fid.write(line)




    # ----------------- 控制参数变量定义 -----------------
    # Figure 尺寸
    FIG3_SIZE = (8, 6)
    FIG5_SIZE = (6, 6)
    FIG4_SIZE = (6, 6)

    # 通用样式
    MARKER_SIZE = 4              # 点图标记大小（Figure 3）
    TICK_LABEL_SIZE = 12        # 坐标轴刻度标签大小（Figure 3、Figure 4）
    TICK_LABEL_SIZE_FIG5 = 12    # 坐标轴刻度标签大小（Figure 5）
    AXIS_LABEL_SIZE = 12         # 坐标轴标签字体大小（Figure 4）
    TITLE_FONT_SIZE = 18         # 标题字体大小（Figure 4）
    LEGEND_FONT_SIZE = 14        # 图例字体大小（Figure 4）
    nFontsize = 16

    # Figure 4 专用
    LINE_WIDTH_CUMULATIVE = 1.5  # 累积曲线线宽（Figure 4）
    MARKER_SIZE_CUMULATIVE = 2   # 累积曲线标记大小（Figure 4）

    # Figure 5 专用
    BIN_WIDTH = 0.1              # 二维直方图的 bin 宽度
    LINE_WIDTH_LINE = 1          # 直线宽度（Figure 5）

    # ----------------- 绘制图形代码 -----------------

    # ----- Figure 3: 原始数据与去相关目录的纬度随时间变化 -----
    fig3, (ax1, ax2) = plt.subplots(2, 1, figsize=FIG3_SIZE)
    ax1.plot(t, data[:, 6], '.', markersize=MARKER_SIZE)
    ax1.set_ylabel('Latitude',fontsize=nFontsize)
    ax1.xaxis_date()
    fig3.autofmt_xdate()
    ax1.tick_params(labelsize=TICK_LABEL_SIZE)

    # 根据去相关目录前 6 列重构时间
    years_out   = Declustered_Catalog[:, 0].astype(int)
    months_out  = Declustered_Catalog[:, 1].astype(int)
    days_out    = Declustered_Catalog[:, 2].astype(int)
    hours_out   = Declustered_Catalog[:, 3].astype(int)
    minutes_out = Declustered_Catalog[:, 4].astype(int)
    secs_out    = Declustered_Catalog[:, 5]
    t1_list = []
    for i in range(len(years_out)):
        t1_list.append(datetime(years_out[i], months_out[i], days_out[i],
                                hours_out[i], minutes_out[i], int(secs_out[i])))
    t1 = mdates.date2num(t1_list)
    ax2.plot(t1, Declustered_Catalog[:, 6], '.', markersize=MARKER_SIZE)
    ax2.set_ylabel('Latitude',fontsize=nFontsize)
    ax2.set_xlabel('Time', fontsize=nFontsize)
    ax2.xaxis_date()
    fig3.autofmt_xdate()
    ax2.tick_params(labelsize=TICK_LABEL_SIZE)
    fig3.tight_layout()
    fig3.savefig(config.decluster_data_fig_lat, dpi=300)

    # ----- Figure 5: 二维直方图及平滑图像 -----
    num_bins = int(round(10 / BIN_WIDTH))
    X = Dataobj[:, 0]
    Y = Dataobj[:, 1]
    A = np.zeros((num_bins, num_bins))
    # 计算每个网格中事件数
    for i in range(num_bins):
        for j in range(num_bins):
            x_min = -9 - BIN_WIDTH + BIN_WIDTH * (i + 1)
            x_max = -9 + BIN_WIDTH * (i + 1)
            y_min = -6 - BIN_WIDTH + BIN_WIDTH * (j + 1)
            y_max = -6 + BIN_WIDTH * (j + 1)
            A[j, i] = np.sum((X >= x_min) & (X < x_max) & (Y > y_min) & (Y < y_max))
    x_extent = [-9, -9 + BIN_WIDTH * num_bins]
    y_extent = [-6, -6 + BIN_WIDTH * num_bins]

    fig5, ax = plt.subplots(figsize=FIG5_SIZE)
    im = ax.imshow(A, extent=(x_extent[0], x_extent[1], y_extent[0], y_extent[1]),
                   origin='lower', cmap='jet')
    plt.colorbar(im, ax=ax)
    # 平滑处理
    Nr, Nc = 4, 4
    matrixOut = mdnna.smooth2a(A, Nr, Nc)
    ax.imshow(matrixOut, extent=(x_extent[0], x_extent[1], y_extent[0], y_extent[1]),
              origin='lower', cmap='jet')
    ax.set_xlim([-7, -1])
    ax.set_ylim([-5, 3])
    # 绘制绿色直线：T + R = constant，即 R = eta0 - T
    T_line = np.arange(-10, 3, 1)
    R_line = eta0 - T_line
    ax.plot(T_line, R_line, color='white', linewidth=LINE_WIDTH_LINE)
    ax.set_xlabel('Rescaled time, T',fontsize=nFontsize)
    ax.set_ylabel('Rescaled distance, R',fontsize=nFontsize)
    ax.set_xticks([-7, -6, -5, -4, -3, -2, -1])
    ax.set_yticks([-5, -4, -3, -2, -1, 0, 1, 2, 3])
    ax.set_xticklabels(
        [r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticklabels(
        [r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$', r'$10^{2}$',
         r'$10^{3}$'])
    ax.tick_params(labelsize=TICK_LABEL_SIZE_FIG5)
    # fig5.tight_layout()
    fig5.savefig(config.decluster_data_fig_nna, dpi=300)

    # ----- Figure 4: 累积地震次数随时间变化 -----
    # 原始目录
    t_orig_list = []
    len_data = data.shape[0]
    for i in range(len_data):
        t_orig_list.append(datetime(int(data[i, 0]), int(data[i, 1]), int(data[i, 2]),
                                    int(data[i, 3]), int(data[i, 4]), int(data[i, 5])))
    t_orig = mdates.date2num(t_orig_list)
    num_eq = np.arange(1, len(t_orig) + 1)
    num_eq_normalized = num_eq / num_eq[-1]
    fig4, ax4 = plt.subplots(figsize=FIG4_SIZE)
    ax4.plot(t_orig, num_eq_normalized, '-o', linewidth=LINE_WIDTH_CUMULATIVE, markersize=MARKER_SIZE_CUMULATIVE,
             color='black', label='Original Catalog')
    # 去相关目录
    t_decl_list = []
    for i in range(Declustered_Catalog.shape[0]):
        t_decl_list.append(datetime(int(Declustered_Catalog[i, 0]), int(Declustered_Catalog[i, 1]),
                                    int(Declustered_Catalog[i, 2]), int(Declustered_Catalog[i, 3]),
                                    int(Declustered_Catalog[i, 4]), int(Declustered_Catalog[i, 5])))
    t_decl = mdates.date2num(t_decl_list)
    num_eq_decl = np.arange(1, len(t_decl) + 1)
    num_eq_decl_normalized = num_eq_decl / num_eq_decl[-1] if len(num_eq_decl) > 0 else num_eq_decl
    ax4.plot(t_decl, num_eq_decl_normalized, '-o', linewidth=LINE_WIDTH_CUMULATIVE, markersize=MARKER_SIZE_CUMULATIVE,
             color='red', label='Declustered Catalog')
    ax4.set_xlabel('Time', fontsize=nFontsize)
    ax4.set_ylabel('Normalized Cumulative Number of Earthquakes', fontsize=nFontsize)
    # ax4.set_title('Normalized Cumulative Earthquake Count Over Time', fontsize=TITLE_FONT_SIZE, fontweight='bold')
    ax4.grid(True)
    ax4.legend(fontsize=LEGEND_FONT_SIZE, loc='upper left')
    ax4.xaxis_date()
    fig4.autofmt_xdate()
    ax4.tick_params(labelsize=TICK_LABEL_SIZE)
    fig4.tight_layout()
    fig4.savefig(config.decluster_data_fig_cum, dpi=300)