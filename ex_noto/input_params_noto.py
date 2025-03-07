import os
from datetime import datetime
import sys

#
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'config'))
print("src_path =", src_path)  # 检查路径是否正确
sys.path.insert(0, src_path)
from config import Config



class ExtendedConfig(Config):
    def __init__(self):

        BASE_PATH = '../ex_noto/'
        Study_region = 'Noto_'
        super().__init__(BASE_PATH,Study_region)

        self.data_original = 'input/JMA_noto_eq.txt'
        self.input_stain = 'input/solid-ocean.out'

        self.mainshock_t = datetime(2023,12,31,0,0,0)

        # 断层参数
        self.str = 213
        self.dip = 41
        self.slip = 79
        self.depth = 6  # km
        self.main_lat = 37.5
        self.main_lon = 137

        #%% 选目录参数
        self.lat_center = 37.5
        self.lon_center = 137.00
        self.side_length = 1
        self.side_length_x = 1.8
        self.side_length_y = 1.2
        # self.start_time = datetime(2003, 1, 1, 0, 0, 0)
        self.start_time = datetime(2006, 1, 1, 0, 0, 0)

        self.end_time = datetime(2023, 12, 31, 0, 0, 0)
        # self.end_time = datetime(2020, 12, 1, 0, 0, 0)

        self.depth_cut = 20
        self.angle_deg = 0
        self.mag_cut = None


        #%% 去相关参数
        self.df = 1.6
        self.b = 1.05
        self.p = 0.5
        self.q = 1 - self.p
        self.eta0 = -4.5

        # 相位分箱: -180 到 180，步长 10
        self.bin_phase = list(range(-180, 180,20))

        #%% 潮汐和应力参数
        self.t_stress_start = '2002-06-01'
        self.t_sample = 360
        self.miu = 0.6
        self.t_search = 3
        self.ref_interval = 1 / 24

        #%% 时间评价参数
        self.window_days = 3 * 365
        self.step_days = 0.1*self.window_days
        self.syn_times = 100
        self.initial_guess = [2, 0.01]