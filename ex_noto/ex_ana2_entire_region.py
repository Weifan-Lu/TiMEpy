# main.py
import sys
import os
# from input_params_rc import config
from input_params_noto import ExtendedConfig

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, src_path)

# import src  # At this point, the src package can be successfully imported

# Analysis
from analysis.ana2_entire_region import ana2_entire_region

def main():
    config = ExtendedConfig()
    ana2_entire_region(config)

if __name__ == '__main__':
    main()
