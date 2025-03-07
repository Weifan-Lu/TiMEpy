import os
from datetime import datetime
import sys


src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'config'))
sys.path.insert(0, src_path)
from config import Config

class ExtendedConfig(Config):
    def __init__(self):
        BASE_PATH = '../ex_noto/'
        Study_region = 'Noto_'
        super().__init__(BASE_PATH, Study_region)

        # Base path
        self.data_original = 'input/JMA_noto_eq.txt'
        self.input_stain = 'input/solid-ocean.out'

        # Mainshock information
        self.mainshock_t = datetime(2023, 12, 31, 0, 0, 0)
        self.str, self.dip, self.slip = 213, 41, 79
        self.main_lat, self.main_lon, self.depth = 37.5, 137.0, 6

        # %% Selection region parameters
        self.lat_center, self.lon_center = 37.5, 137.00
        self.side_length_x, self.side_length_y, self.angle_deg = 1.8, 1.2, 0
        self.start_time, self.end_time = datetime(2006, 1, 1, 0, 0, 0), datetime(2023, 12, 31, 0, 0, 0)
        self.depth_cut = 20
        self.mag_cut = None

        # %% Decorrelation parameters
        self.df, self.b, self.p = 1.6, 1.05, 0.5
        self.q = 1 - self.p
        self.eta0 = -4.5

        # Phase binning: -180 to 180, step size 20
        self.bin_phase = list(range(-180, 181, 20))

        # %% Tidal and stress parameters
        self.t_stress_start = '2002-06-01'
        self.t_sample = 360
        self.miu = 0.6
        self.t_search = 3
        self.ref_interval = 1 / 24

        # %% Temporal evaluation parameters
        self.window_days, self.step_days = 3 * 365, 0.1 * (3 * 365)
        self.initial_guess = [2, 0.01]
