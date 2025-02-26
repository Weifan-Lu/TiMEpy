import os
from datetime import datetime
import sys

#
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'config'))
print("src_path =", src_path)
sys.path.insert(0, src_path)
from config import Config



class ExtendedConfig(Config):
    def __init__(self):
        BASE_PATH = '../ex_ridgecrest/'
        Study_region = 'New_RC_'
        super().__init__(BASE_PATH,Study_region)
        # Base path
        self.data_original = 'input/ridgrest_10years.txt'

        # Mainshock time
        self.mainshock_t = datetime(2019, 7, 4, 0, 0, 0)

        # Fault parameters
        self.str = 321
        self.dip = 81
        self.slip = 180
        self.depth = 12  # km
        self.main_lat = 35.76950
        self.main_lon = -117.59933

        # %% Selection region parameters
        self.lat_center = 35.72
        self.lon_center = -117.56
        self.side_length = 0.4
        self.side_length_x = 0.1
        self.side_length_y = 1.0
        self.start_time = datetime(2010, 1, 1, 0, 0, 0)
        self.end_time = datetime(2019, 7, 3, 0, 0, 0)
        self.depth_cut = 30
        self.angle_deg = 135

        # %% Decorrelation parameters
        self.df = 1.6
        self.b = 1.05
        self.p = 0.5
        self.q = 1 - self.p
        self.eta0 = -4.5

        # Phase binning: -180 to 180, step size 10
        self.bin_phase = list(range(-180, 181, 20))

        # %% Tidal and stress parameters
        self.t_stress_start = '2009-01-01'
        self.t_sample = 360
        self.miu = 0.6
        self.t_search = 3
        self.ref_interval = 1 / 24

        # %% Temporal evaluation parameters
        self.window_days = 3.4 * 365
        self.step_days = 91
        self.syn_times = 100
        self.initial_guess = [2, 0.1]
