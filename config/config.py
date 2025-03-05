import os

class Config:
    def __init__(self,BASE_PATH,Study_region):

        self.output = os.path.join(BASE_PATH, 'output')
        self.output_catalog = os.path.join(BASE_PATH, 'output', 'catalog')
        self.output_tidal_phase = os.path.join(BASE_PATH, 'output', 'tidal_phase')
        self.output_stress = os.path.join(BASE_PATH, 'output', 'stress')
        self.output_fig = os.path.join(BASE_PATH, 'output', 'figure')

        self.input = os.path.join(BASE_PATH, 'input')
        self.data_select = os.path.join(self.output_catalog, Study_region + 'select_catalog.txt')

        self.data_select_decluster = os.path.join(self.output_catalog, Study_region + 'select_catalog_decluster.txt')
        self.data_select_decluster_dataobj = os.path.join(self.output_catalog, Study_region + 'select_catalog_decluster_dataobj.npy')

        self.phase_stress_obs = os.path.join(self.output_tidal_phase, Study_region + 'phase_stress_obs')
        self.phase_stress_ref = os.path.join(self.output_tidal_phase, Study_region + 'phase_stress_ref')

        self.input_stain = os.path.join(self.input, 'solid-ocean.out')
        self.output_stress_txt = os.path.join(self.output_stress,Study_region + 'stress_Vol_N_S.txt')

        # 图形文件
        self.select_data_fig = os.path.join(self.output_catalog, Study_region + 'select_catalog.png')
        self.select_data_fig_map = os.path.join(self.output_catalog, Study_region + 'select_catalog_map.pdf')
        self.decluster_data_fig_lat = os.path.join(self.output_catalog, Study_region + 'decluster_catalog_lat.png')
        self.decluster_data_fig_cum = os.path.join(self.output_catalog, Study_region + 'decluster_catalog_cum.png')
        self.decluster_data_fig_nna = os.path.join(self.output_catalog, Study_region + 'decluster_catalog_nna.png')

        self.TM_all_region_Vol = os.path.join(self.output_fig, Study_region + 'entire_region_Vol.pdf')
        self.TM_all_region_S = os.path.join(self.output_fig, Study_region+ 'entire_region_S.pdf')
        self.TM_all_region_N = os.path.join(self.output_fig, Study_region + 'entire_region_N.pdf')
        self.TM_all_region_CFS = os.path.join(self.output_fig, Study_region + 'entire_region_CFS.pdf')
        self.TM_all_region_tidal_sensitivity = os.path.join(self.output_fig, Study_region + 'TM_all_region_tidal_sensitivity.pdf')

        self.TM_all_region_tidal_sensitivity_vol = os.path.join(self.output_fig, Study_region + 'TM_all_region_tidal_sensitivity_vol.pdf')
        self.TM_all_region_tidal_sensitivity_n = os.path.join(self.output_fig, Study_region + 'TM_all_region_tidal_sensitivity_n.pdf')
        self.TM_all_region_tidal_sensitivity_s = os.path.join(self.output_fig, Study_region + 'TM_all_region_tidal_sensitivity_s.pdf')
        self.TM_all_region_tidal_sensitivity_cfs = os.path.join(self.output_fig, Study_region + 'TM_all_region_tidal_sensitivity_cfs.pdf')

        self.TM_time_tidal_P_value = os.path.join(self.output_fig, Study_region + 'temporal_variation_P_value.pdf')
        self.TM_time_tidal_sensitivity = os.path.join(self.output_fig, Study_region + 'temporal_variation_tidal_sensitivity.pdf')
        self.TM_time_tidal_phase_shift = os.path.join(self.output_fig, Study_region + 'temporal_variation_tidal_phase_shift.pdf')
        self.TM_time_tidal_amplitude = os.path.join(self.output_fig, Study_region+ 'temporal_variation_tidal_amplitude.pdf')
        self.TM_time_earthquake_number = os.path.join(self.output_fig, Study_region + 'temporal_variation_earthquake_number.pdf')

        self.TM_all_region_stress = os.path.join(self.output_fig, Study_region + 'TM_time_tidal_stress.pdf')
        self.TM_time_tidal_all_CFS = os.path.join(self.output_fig, Study_region + 'temporal_variation_CFS.pdf')
        self.TM_earthquake_b_positive = os.path.join(self.output_fig, Study_region + 'b_value.pdf')
