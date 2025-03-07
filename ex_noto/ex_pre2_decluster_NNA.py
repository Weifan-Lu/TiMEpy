# main.py
import sys
import os
# from input_params_rc import config
from input_params_noto import ExtendedConfig

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, src_path)

# import src  # At this point, the src package can be successfully imported

from preprocessing.pre2_decluster_NNA import pre2_decluster_NNA


def main():
    config = ExtendedConfig()
    pre2_decluster_NNA(config)


if __name__ == '__main__':
    main()
