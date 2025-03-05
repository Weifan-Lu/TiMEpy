import numpy as np
import math
from datetime import datetime, timedelta
import main_strain_to_stress as sts


# --- Main Function: Convert Strain to Stress ---
def pre3_strain_to_stress(config):
    """
    Reads relevant files and parameters from the parameter dictionary `config`,
    calculates stress changes, and saves the results as a .mat file and a TXT file.

    The `config` dictionary must include:
      - input_stain: Path to the strain data file
      - output_stress_txt: Path to the output stress result text file
      - miu: Friction coefficient (fc)
      - str: Fault strike angle (degrees)
      - dip: Fault dip angle (degrees)
      - slip: Fault slip angle (degrees)
      - depth: Depth (km)
      - stress_to_mat: Path to save the stress result in a .mat file
      - t_stress_start: Start time for stress calculation (datetime object or string in 'YYYY-MM-DD' format)
      - t_sample: Sampling interval (seconds)
    """

    print('====== Processing | Strin to Stress: Start ======')

    infile = config.input_stain
    outfile = config.output_stress_txt
    fc = config.miu
    str_val = config.str
    dip = config.dip
    slip = config.slip
    depth = config.depth

    t_stress_start = config.t_stress_start  # Should be a datetime object or a string in 'YYYY-MM-DD' format
    t_sample = config.t_sample

    # Compute Lamé parameters
    alamda, amiu, depth, vp, vs, rho = sts.select_Lame(depth)

    # Read strain data and compute stress changes
    normal_stress, shear_stress, volumetric_strain, t_tide = sts.strain_to_stress(
        infile, alamda, amiu, fc, str_val, dip, slip, t_stress_start, t_sample
    )

    # Convert datetime objects in t_tide to string format before saving
    t_tide_str = np.array([dt.strftime('%Y-%m-%d %H:%M:%S') for dt in t_tide])
    # savemat(stress_to_mat, {
    #     'normal_stress': normal_stress,
    #     'shear_stress': shear_stress,
    #     'volumetric_strain': volumetric_strain,
    #     't_tide': t_tide_str
    # })

    # Read original strain file (assuming the first two columns contain date and time, and the 9th column contains a numerical value)
    with open(infile, 'r') as f:
        orig_lines = f.readlines()

    # Write results to the output text file
    with open(outfile, 'w') as f_out:
        neq = len(normal_stress)
        for i in range(neq):
            tokens = orig_lines[i].strip().split()
            # Skip if the line does not have enough information
            if len(tokens) < 9:
                continue
            # Extract original date and time
            date_str = tokens[0]
            time_str = tokens[1]
            # Convert to string if they are datetime objects (usually, the original data is already in string format)
            if isinstance(date_str, datetime):
                date_str = date_str.strftime('%Y-%m-%d')
            if isinstance(time_str, datetime):
                time_str = time_str.strftime('%H:%M:%S')
            try:
                col9 = float(tokens[8])
            except Exception as e:
                col9 = 0.0
            # Write output in the format: Date Time  Original 9th column value  Volumetric Strain  Shear Stress  Normal Stress
            f_out.write(
                f"{date_str} {time_str}  {col9:.6f}  {volumetric_strain[i]:.8e}  {shear_stress[i]:.8e}  {normal_stress[i]:.8e}\n"
            )

    print(f"Stress results have been written to {outfile}")
    print('====== Processing | Strin to Stress: End ======')
