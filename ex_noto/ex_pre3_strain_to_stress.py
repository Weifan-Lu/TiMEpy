# main.py
import sys
import os
# from input_params_rc import config
from input_params_noto import ExtendedConfig

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, src_path)

# import src  # At this point, the src package can be successfully imported

from preprocessing.pre3_strain_to_stress import pre3_strain_to_stress

def main():
    config = ExtendedConfig()
    pre3_strain_to_stress(config)

if __name__ == '__main__':
    main()
