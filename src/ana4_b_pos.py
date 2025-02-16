import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from b_posivite import b_positive_discrete, b_positive_discrete_bootstrap # 请确保该模块和函数能正常调用
import matplotlib as mpl

# 设置全局字体大小为12
mpl.rcParams['font.size'] = 12
def ana4_b_pos(config):
    plt.close('all')

    # 1. 加载数据
    data = np.loadtxt(config.data_select)

    # 2. 按时间排序（假设数据的前6列依次为 年, 月, 日, 时, 分, 秒）
    sort_indices = np.lexsort((data[:, 5], data[:, 4], data[:, 3],
                               data[:, 2], data[:, 1], data[:, 0]))
    sorted_data = data[sort_indices]

    # 3. 提取震级数据（假定震级位于第10列，即索引9）
    magnitudes = sorted_data[:, 9]

    # 4. 将时间数据转换为 datetime 对象（这里假设前6列依次为 年, 月, 日, 时, 分, 秒）
    t_data = np.array([
        datetime.datetime(int(row[0]), int(row[1]), int(row[2]),
                          int(row[3]), int(row[4]), int(row[5]))
        for row in sorted_data
    ])

    # 5. 获取时间范围和滑动窗口参数
    dt_start = config.start_time  # 应为 datetime 对象
    dt_end = config.end_time  # 应为 datetime 对象
    tmin = mdates.date2num(dt_start)
    tmax = mdates.date2num(dt_end)
    step_days = config.step_days  # 步长（单位：天）
    window_days = config.window_days  # 窗口长度（单位：天）

    # 准备存储结果的列表
    window_centers = []
    b_values = []
    b_std_all = []
    # 6. 滑动窗口遍历
    # np.arange 生成的 t_start 为 matplotlib 的日期数字
    t_data_num = mdates.date2num(t_data)  # 便于比较
    for t_start in np.arange(tmin, tmax - window_days + 1e-6, step_days):
        t_end = t_start + window_days
        # 选择当前窗口内的观测数据
        indices = np.where((t_data_num >= t_start) & (t_data_num <= t_end))[0]
        mag_cut = magnitudes[indices]

        # 若窗口内有足够数据，则计算 b 值，否则设为 NaN
        if len(mag_cut) > 0:
            b_val = b_positive_discrete(mag_cut, delta=0.05)
            b_mean, b_std = b_positive_discrete_bootstrap(magnitudes, delta=0.05, n_boot=1000)
        else:
            b_val = np.nan

        # 保存窗口中心的时间（这里采用窗口起始时间和结束时间的中点）
        window_center = t_start + window_days
        window_centers.append(window_center)
        b_values.append(b_val)
        b_std_all.append(b_std)

    # 将窗口中心从数字格式转换为 datetime 格式
    window_dates = mdates.num2date(window_centers)

    # 创建两个子图，共享 x 轴
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(14, 6))
    # 子图1：绘制带误差棒的 b 值随时间变化图
    ax1.errorbar(
        window_dates, b_values, yerr=b_std_all, fmt='o',
        markersize=6, markerfacecolor='royalblue', markeredgecolor='black',  # 美化标记
        capsize=4, capthick=1.2, elinewidth=1, alpha=0.8, color='royalblue',  # 美化误差棒
        linestyle='-', label='b-value'
    )
    ax1.set_ylabel('b-value')  # 字体大小已全局设置
    # ax1.set_title('Time vs. b-value')
    ax1.legend()

    # 子图2：绘制地震震级随时间变化图，震级用灰色圆圈点表示
    ax2.plot(
        t_data, magnitudes, 'o',
        markersize=4, color='gray', alpha=0.8, label='Magnitude'
    )
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Magnitude')
    # ax2.set_title('Time vs. Magnitude')
    ax2.legend()

    # 设置 x 轴日期格式（共享 x 轴）
    ax2.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax2.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax2.xaxis.get_major_locator()))
    fig.autofmt_xdate()

    plt.tight_layout()
    plt.savefig(config.TM_earthquake_b_positive, format='pdf')
    plt.show()