import numpy as np
import math
from datetime import datetime, timedelta

def emprem(rr, ind):
    """
    Calculate physical parameters at a given radius based on the PREM Earth model.

    Parameters:
      rr  : Radius (km), 0 ≤ rr ≤ 6371
      ind : Adjustment factor, if ind < 0 and the layer is not the outermost, decrement the layer index by 1

    Returns:
      rho : Density
      vp  : P-wave velocity
      vs  : S-wave velocity
      qp  : P-wave quality factor
      qs  : S-wave quality factor
    """
    # Define PREM model parameters
    r = np.array([0, 1221.5, 3480, 3630, 5600, 5701, 5771, 5971, 6151, 6291, 6346.6, 6356, 6368, 6371], dtype=float)
    # d: 4 x 13 array, each column corresponds to parameters of a layer
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

    # Check if rr is within the valid range
    if rr < 0 or rr > 6371:
        raise ValueError("Radius is outside the valid range for the PREM model.")

    # Find the layer: locate the last j satisfying r[j] <= rr < r[j+1]
    i = None
    for j in range(len(r) - 1):
        if rr >= r[j] and rr < r[j + 1]:
            i = j
    if i is None and np.isclose(rr, r[-1]):
        i = len(r) - 2  # Select the outermost layer
    if i is None:
        raise ValueError("Radius does not fall within any layer in the PREM model.")

    # If ind < 0 and i > 0, adjust the layer index
    if ind < 0 and i > 0:
        i = i - 1

    # Normalize radius
    x = rr / r[-1]
    print(f"rr = {rr}, i = {i}")

    # Calculate density, P-wave velocity, and S-wave velocity (using Horner's method)
    rho = d[0, i] + x * (d[1, i] + x * (d[2, i] + x * d[3, i]))
    vp = p[0, i] + x * (p[1, i] + x * (p[2, i] + x * p[3, i]))
    vs = s[0, i] + x * (s[1, i] + x * (s[2, i] + x * s[3, i]))

    # S-wave quality factor
    qs_val = qm[i]
    if qs_val > 1e10 or np.isnan(qs_val) or np.isinf(qs_val):
        qs_val = 0

    # Calculate P-wave quality factor
    vsvp = vs / vp
    vsvp2 = vsvp ** 2
    al = vsvp2 * 4 / 3
    qpinv = al / qm[i] + (1 - al) / qk[i]
    qp = 1 / qpinv

    return rho, vp, vs, qp, qs_val


def emdlv(r_val, ind):
    """
    Earth model setup for PREM
    Directly calls the emprem function to calculate parameters.

    Parameters:
      r_val : Radius (km)
      ind   : Adjustment factor
    Returns:
      Same five parameters as emprem
    """
    return emprem(r_val, ind)


# --- Calculate Lamé parameters ---
def select_Lame(depth):
    """
    Calculate Lamé parameters (λ and μ) and other parameters based on the given depth
    """
    r = 6371.0 - depth  # Radius in km
    ind = 0
    rho, vp, vs, qp, qs = emdlv(r, ind)
    lam = rho * (vp ** 2 - 2 * vs ** 2)
    mu = rho * vs ** 2
    print(f"{lam:.1f} {mu:.1f} {depth:.3f} {vp:.2f} {vs:.2f} {rho:.2f}")
    return lam, mu, depth, vp, vs, rho


# --- Calculate stress changes from strain tensor ---
def strn2stress(strn, alamda, amiu, fc, str_val, dip, slip):
    """
    Input:
      strn: 6 x n numpy array, representing the 6 components of the strain tensor (row order: exx, eyy, ezz, exy, exz, eyz)
      alamda, amiu: Lamé parameters (λ and μ)
      fc: Friction coefficient
      str_val, dip, slip: Fault strike, dip, and slip angles (in degrees)
    Output:
      dT: Shear stress change (sgmxy)
      dS: Normal stress change (sgmyy)
      dCFF: Coulomb stress change = sgmxy + fc * sgmyy
    """
    n = strn.shape[1]
    dT = np.zeros(n)
    dS = np.zeros(n)
    dCFF = np.zeros(n)

    # Convert angles to radians
    rad = math.pi / 180.0
    strr = str_val * rad
    dipr = dip * rad
    slipr = slip * rad

    # Calculate trigonometric values
    css = math.cos(strr)
    sns = math.sin(strr)
    csd = math.cos(dipr)
    snd = math.sin(dipr)
    csl = math.cos(slipr)
    snl = math.sin(slipr)

    # Calculate fault slip vector and normal vector components
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


# --- Calculate stress changes from strain file ---
def strain_to_stress(infile, alamda, amiu, fc, str_val, dip, slip, t_stress_start, t_sample):
    """
    Read strain data file, calculate volumetric strain and stress changes.
    Each line in the file contains 7 floating-point numbers starting from the 20th character,
    where the first 6 are strain tensor components, and the 7th can be ignored.
    Volumetric strain dV = exx + eyy + ezz.

    Returns:
      normal_stress: Normal stress change (dS)
      shear_stress: Shear stress change (dT)
      volumetric_strain: Volumetric strain change dV
      t_tide: Time vector starting from t_stress_start with sampling interval t_sample (seconds) (numpy array, each element is a datetime object)
    """
    lines = []
    with open(infile, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                lines.append(line)
    neq = len(lines)
    # Preallocate strain tensor and volumetric strain arrays
    strn = np.zeros((6, neq))
    dV = np.zeros(neq)

    for i, line in enumerate(lines):
        parts = line.split()  # Split by spaces to get columns
        try:
            # Assume data columns start from the third column, take 6 values (strain data)
            # parts[0] is date, parts[1] is time, parts[2:8] are strain data
            numbers = list(map(float, parts[2:8]))
            if len(numbers) < 6:
                raise ValueError("Insufficient data")
        except Exception as e:
            raise ValueError(f"Error parsing line {i + 1}: {line}\n{e}")

        strn[:, i] = numbers
        # Calculate volumetric strain: sum of the first three components
        dV[i] = numbers[0] + numbers[1] + numbers[2]

    # Call strn2stress to calculate stress changes
    dT, dS, dCFF = strn2stress(strn, alamda, amiu, fc, str_val, dip, slip)

    # Assume t_stress_start is a string
    if isinstance(t_stress_start, str):
        # print(t_stress_start)
        t_stress_start = datetime.strptime(t_stress_start, '%Y-%m-%d')

    # Generate time vector t_tide, starting from t_stress_start, with sampling interval t_sample seconds
    t_tide = np.array([t_stress_start + timedelta(seconds=int(i * t_sample)) for i in range(neq)])

    normal_stress = dS * 1e9  # Normal stress change
    shear_stress = dT * 1e9  # Shear stress change
    volumetric_strain = dV * 1e10

    return normal_stress, shear_stress, volumetric_strain, t_tide