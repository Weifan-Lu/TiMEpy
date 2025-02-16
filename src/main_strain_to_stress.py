import numpy as np
import math
from datetime import datetime, timedelta


def emprem(rr, ind):
    """
    根据 PREM 地球模型计算给定半径下的物理参数。

    参数:
      rr  : 半径（km），0 ≤ rr ≤ 6371
      ind : 调整因子，如果 ind < 0 并且所在层不在最外层，则将层编号减 1

    返回:
      rho : 密度
      vp  : P 波速度
      vs  : S 波速度
      qp  : P 波品质因子
      qs  : S 波品质因子
    """
    # 定义 PREM 模型参数
    r = np.array([0, 1221.5, 3480, 3630, 5600, 5701, 5771, 5971, 6151, 6291, 6346.6, 6356, 6368, 6371], dtype=float)
    # d: 4 x 13 数组，每列对应一个层的参数
    d = np.array([[13.0885, 12.5815, 7.9565, 7.9565, 7.9565, 5.3197, 11.2494, 7.1089, 2.691, 2.691, 2.9, 2.6, 2.6],
                  [0, -1.2638, -6.4761, -6.4761, -6.4761, -1.4836, -8.0298, -3.8045, 0.6924, 0.6924, 0, 0, 0],
                  [-8.8381, -3.6426, 5.5283, 5.5283, 5.5283, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, -5.5281, -3.0807, -3.0807, -3.0807, 0, 0, 0, 0, 0, 0, 0, 0]])
    p = np.array([[11.2622, 11.0487, 15.3891, 24.952, 29.2766, 19.0957, 39.7027, 20.3926, 4.1875, 4.1875, 6.8, 5.8, 5],
                  [0, -4.0362, -5.3181, -40.4673, -23.6027, -9.8672, -32.6166, -12.2569, 3.9382, 3.9382, 0, 0, 0],
                  [-6.364, 4.8023, 5.5242, 51.4832, 5.5242, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, -13.5732, -2.5514, -26.6419, -2.5514, 0, 0, 0, 0, 0, 0, 0, 0]])
    s = np.array([[3.6678, 0, 6.9254, 11.1671, 22.3459, 9.9839, 22.3512, 8.9496, 2.1519, 2.1519, 3.9, 3.2, 2.6],
                  [0, 0, 1.4672, -13.7818, -17.2473, -4.9324, -18.5856, -4.4597, 2.3481, 2.3481, 0, 0, 0],
                  [-4.4475, 0, -2.0834, 17.4575, -2.0834, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0.9783, -9.2777, 0.9783, 0, 0, 0, 0, 0, 0, 0, 0]])
    qm = np.array([84.6, 1e30, 312, 312, 312, 143, 143, 143, 80, 600, 600, 600, 600])
    qk = np.array([1327.7, 57823, 57823, 57823, 57823, 57823, 57823, 57823, 57823, 57823, 57823, 57823, 57823])

    # 检查 rr 是否在有效范围内
    if rr < 0 or rr > 6371:
        raise ValueError("Radius is outside the valid range for the PREM model.")

    # 查找所在层：找到最后一个满足 r[j] <= rr < r[j+1] 的 j
    i = None
    for j in range(len(r) - 1):
        if rr >= r[j] and rr < r[j + 1]:
            i = j
    if i is None and np.isclose(rr, r[-1]):
        i = len(r) - 2  # 选取最外层
    if i is None:
        raise ValueError("Radius does not fall within any layer in the PREM model.")

    # 如果 ind < 0 并且 i > 0，则调整层编号
    if ind < 0 and i > 0:
        i = i - 1

    # 归一化半径
    x = rr / r[-1]
    print(f"rr = {rr}, i = {i}")

    # 计算密度、P 波速度和 S 波速度（使用 Horner 法则）
    rho = d[0, i] + x * (d[1, i] + x * (d[2, i] + x * d[3, i]))
    vp = p[0, i] + x * (p[1, i] + x * (p[2, i] + x * p[3, i]))
    vs = s[0, i] + x * (s[1, i] + x * (s[2, i] + x * s[3, i]))

    # S 波品质因子
    qs_val = qm[i]
    if qs_val > 1e10 or np.isnan(qs_val) or np.isinf(qs_val):
        qs_val = 0

    # 计算 P 波品质因子
    vsvp = vs / vp
    vsvp2 = vsvp ** 2
    al = vsvp2 * 4 / 3
    qpinv = al / qm[i] + (1 - al) / qk[i]
    qp = 1 / qpinv

    return rho, vp, vs, qp, qs_val


def emdlv(r_val, ind):
    """
    Earth model setup for PREM
    直接调用 emprem 函数计算参数。

    参数:
      r_val : 半径（km）
      ind   : 调整因子
    返回:
      与 emprem 相同的五个参数
    """
    return emprem(r_val, ind)


# --- 计算 Lamé 参数 ---
def select_Lame(depth):
    """
    根据给定深度计算 Lamé 参数（λ 和 μ）及其它参数
    """
    r = 6371.0 - depth  # Radius in km
    ind = 0
    rho, vp, vs, qp, qs = emdlv(r, ind)
    lam = rho * (vp ** 2 - 2 * vs ** 2)
    mu = rho * vs ** 2
    print(f"{lam:.1f} {mu:.1f} {depth:.3f} {vp:.2f} {vs:.2f} {rho:.2f}")
    return lam, mu, depth, vp, vs, rho


# --- 根据应变张量计算应力变化 ---
def strn2stress(strn, alamda, amiu, fc, str_val, dip, slip):
    """
    输入:
      strn: 6 x n numpy 数组，表示应变张量的 6 个分量（行顺序：exx, eyy, ezz, exy, exz, eyz）
      alamda, amiu: Lamé 参数（λ 和 μ）
      fc: 摩擦系数
      str_val, dip, slip: 分别为断层走向、倾角、滑动角（单位：度）
    输出:
      dT: 剪应力变化（sgmxy）
      dS: 正应力变化（sgmyy）
      dCFF: Coulomb 应力变化 = sgmxy + fc * sgmyy
    """
    n = strn.shape[1]
    dT = np.zeros(n)
    dS = np.zeros(n)
    dCFF = np.zeros(n)

    # 角度转换为弧度
    rad = math.pi / 180.0
    strr = str_val * rad
    dipr = dip * rad
    slipr = slip * rad

    # 计算三角函数值
    css = math.cos(strr)
    sns = math.sin(strr)
    csd = math.cos(dipr)
    snd = math.sin(dipr)
    csl = math.cos(slipr)
    snl = math.sin(slipr)

    # 计算断层滑动向量和法向量分量
    a11 = csl * sns - csd * snl * css
    a12 = csl * css + csd * snl * sns
    a13 = snl * snd
    a21 = snd * css
    a22 = -snd * sns
    a23 = csd

    for i in range(n):
        exx = strn[0, i]
        eyy = strn[1, i]
        ezz = strn[2, i]
        exy = strn[3, i]
        exz = strn[4, i]
        eyz = strn[5, i]
        evol = exx + eyy + ezz

        b1 = a21 * exx + a22 * exy + a23 * exz
        b2 = a21 * exy + a22 * eyy + a23 * eyz
        b3 = a21 * exz + a22 * eyz + a23 * ezz

        exy2 = a11 * b1 + a12 * b2 + a13 * b3
        eyy2 = a21 * b1 + a22 * b2 + a23 * b3

        sgmxy = 2 * amiu * exy2
        sgmyy = alamda * evol + 2 * amiu * eyy2

        dT[i] = sgmxy
        dS[i] = sgmyy
        dCFF[i] = sgmxy + fc * sgmyy

    return dT, dS, dCFF


# --- 从应变文件计算应力变化 ---
def strain_to_stress(infile, alamda, amiu, fc, str_val, dip, slip, t_stress_start, t_sample):
    """
    读取应变数据文件，计算体积应变和应力变化。
    文件中每一行从第20个字符开始包含 7 个浮点数，其中前6个为应变张量分量，
    第7个数可忽略。体积应变 dV = exx+eyy+ezz。

    返回:
      normal_stress: 正应力变化（dS）
      shear_stress: 剪应力变化（dT）
      volumetric_strain: 体积应变变化 dV
      t_tide: 以 t_stress_start 为起始、采样间隔为 t_sample（秒）的时间向量（numpy 数组，每个元素为 datetime 对象）
    """
    lines = []
    with open(infile, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                lines.append(line)
    neq = len(lines)
    # 预分配应变张量和体积应变数组
    strn = np.zeros((6, neq))
    dV = np.zeros(neq)

    for i, line in enumerate(lines):
        parts = line.split()  # 按空格分割，得到各列
        try:
            # 假设数据列从第三列开始，取出6个数值（应变数据）
            # parts[0] 为日期, parts[1] 为时间，parts[2:8] 为应变数据
            numbers = list(map(float, parts[2:8]))
            if len(numbers) < 6:
                raise ValueError("数据不足")
        except Exception as e:
            raise ValueError(f"第 {i + 1} 行解析错误: {line}\n{e}")

        strn[:, i] = numbers
        # 计算体积应变：前三个分量之和
        dV[i] = numbers[0] + numbers[1] + numbers[2]

    # 调用 strn2stress 计算应力变化
    dT, dS, dCFF = strn2stress(strn, alamda, amiu, fc, str_val, dip, slip)

    # 假设 t_stress_start 是字符串
    if isinstance(t_stress_start, str):
        # print(t_stress_start)
        t_stress_start = datetime.strptime(t_stress_start, '%Y-%m-%d')

    # 生成时间向量 t_tide，从 t_stress_start 开始，采样间隔为 t_sample 秒
    t_tide = np.array([t_stress_start + timedelta(seconds=int(i * t_sample)) for i in range(neq)])

    normal_stress = dS * 1e9  # 正应力变化
    shear_stress = dT * 1e9  # 剪应力变化
    volumetric_strain = dV * 1e10

    return normal_stress, shear_stress, volumetric_strain, t_tide


