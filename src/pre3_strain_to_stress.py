import numpy as np
import math
from datetime import datetime, timedelta
import main_strain_to_stress as sts


# --- 主函数：应变转换为应力 ---
def pre3_strain_to_stress(config):
    """
    从参数字典 params 中读取相关文件及参数，计算应力变化，并将结果保存为 .mat 文件以及 TXT 文件。

    params 字典中需包含：
      - input_stain: 应变数据文件路径
      - output_stress_txt: 应力结果输出文本路径
      - miu: 摩擦系数（fc）
      - str: 断层走向角（度）
      - dip: 断层倾角（度）
      - slip: 断层滑动角（度）
      - depth: 深度（km）
      - stress_to_mat: 应力结果保存的 .mat 文件路径
      - t_stress_start: 应力计算起始时间（datetime 对象或 'YYYY-MM-DD' 格式字符串）
      - t_sample: 采样间隔（秒）
    """
    infile = config.input_stain
    outfile = config.output_stress_txt
    fc = config.miu
    str_val = config.str
    dip = config.dip
    slip = config.slip
    depth = config.depth

    t_stress_start = config.t_stress_start  # 应为 datetime 对象或 'YYYY-MM-DD' 格式字符串
    t_sample = config.t_sample

    # 计算 Lamé 参数
    alamda, amiu, depth, vp, vs, rho = sts.select_Lame(depth)

    # 读取应变数据并计算应力变化
    normal_stress, shear_stress, volumetric_strain, t_tide = sts.strain_to_stress(infile, alamda, amiu, fc, str_val, dip, slip, t_stress_start, t_sample
    )

    # 在保存之前将 t_tide 中的 datetime 对象转换为字符串
    t_tide_str = np.array([dt.strftime('%Y-%m-%d %H:%M:%S') for dt in t_tide])
    # savemat(stress_to_mat, {
    #     'normal_stress': normal_stress,
    #     'shear_stress': shear_stress,
    #     'volumetric_strain': volumetric_strain,
    #     't_tide': t_tide_str
    # })

    # 读取原始应变文件（假设每行前两列分别为日期和时间，第9列为一个数值）
    with open(infile, 'r') as f:
        orig_lines = f.readlines()

    # 将结果写入输出文本文件
    with open(outfile, 'w') as f_out:
        neq = len(normal_stress)
        for i in range(neq):
            tokens = orig_lines[i].strip().split()
            # 若行内信息不足，则跳过
            if len(tokens) < 9:
                continue
            # 提取原始的日期和时间
            date_str = tokens[0]
            time_str = tokens[1]
            # 如果读取到的是 datetime 对象，则转换为字符串（通常原始数据已为字符串）
            if isinstance(date_str, datetime):
                date_str = date_str.strftime('%Y-%m-%d')
            if isinstance(time_str, datetime):
                time_str = time_str.strftime('%H:%M:%S')
            try:
                col9 = float(tokens[8])
            except Exception as e:
                col9 = 0.0
            # 写入时按照以下格式：日期 时间  原始第9列数值  体积应变  剪应力  正应力
            f_out.write(
                f"{date_str} {time_str}  {col9:.6f}  {volumetric_strain[i]:.8e}  {shear_stress[i]:.8e}  {normal_stress[i]:.8e}\n"
            )

    print(f"应力结果已写入文件 {outfile}")



