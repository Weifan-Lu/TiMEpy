# main.py
import sys
import os
# from input_params_rc import config
from input_params_rc import ExtendedConfig

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
print("src_path =", src_path)
sys.path.insert(0, src_path)

# import src  # At this point, the src package can be successfully imported

# Preprocessing
from preprocessing.pre0_creat_output import pre0_creat_output
from preprocessing.pre1_select_catalog_fault import pre1_select_catalog_fault
from preprocessing.pre2_decluster_NNA import pre2_decluster_NNA
from preprocessing.pre3_strain_to_stress import pre3_strain_to_stress

# Analysis
from analysis.ana1_calc_tidal_phase import ana1_calc_tidal_phase
from analysis.ana2_all_region import ana2_all_region
from analysis.ana3_time_evl import ana3_time_evl
from analysis.ana4_b_pos import ana4_b_pos

def main():
    config = ExtendedConfig()
    print(config)
    # Step 0: Create output directories
    pre0_creat_output(config)
    pre1_select_catalog_fault(config)
    pre2_decluster_NNA(config)
    pre3_strain_to_stress(config)
    opt = '1'
    ana1_calc_tidal_phase(config, opt)
    ana2_all_region(config)
    ana3_time_evl(config)
    ana4_b_pos(config)

if __name__ == '__main__':
    main()
