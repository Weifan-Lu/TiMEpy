# main.py
import sys
import os
# from input_params_rc import config
from input_params_rc import ExtendedConfig

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, src_path)

# import src  # At this point, the src package can be successfully imported

# Preprocessing
from analysis.ana3_temp_variation import ana3_temp_variation

def main():
    config = ExtendedConfig()
    ana3_temp_variation(config)

if __name__ == '__main__':
    main()
