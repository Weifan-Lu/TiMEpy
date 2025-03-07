# main.py
import sys
import os
# from input_params_rc import config
from input_params_noto import ExtendedConfig

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, src_path)


# Analysis
from analysis.ana1_calc_tidal_phase import ana1_calc_tidal_phase

def main():
    config = ExtendedConfig()
    opt = '1'
    ana1_calc_tidal_phase(config, opt)


if __name__ == '__main__':
    main()
