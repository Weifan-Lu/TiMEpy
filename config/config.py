import os

class Config:
    def __init__(self, BASE_PATH, study_region):
        # Base directories
        self.output = os.path.join(BASE_PATH, 'output')
        self.catalog_output = os.path.join(BASE_PATH, 'output', 'catalog')
        self.tidal_phase_output = os.path.join(BASE_PATH, 'output', 'tidal_phase')
        self.stress_output = os.path.join(BASE_PATH, 'output', 'stress')
        self.figure_output = os.path.join(BASE_PATH, 'output', 'figure')

        # Input directory and data files
        self.input = os.path.join(BASE_PATH, 'input')
        self.select_catalog = os.path.join(self.catalog_output, study_region + '_select_catalog.txt')
        self.select_catalog_decluster = os.path.join(self.catalog_output, study_region + '_select_catalog_decluster.txt')
        self.phase_stress_obs = os.path.join(self.tidal_phase_output, study_region + '_phase_stress_obs')
        self.phase_stress_ref = os.path.join(self.tidal_phase_output, study_region + '_phase_stress_ref')

        self.input_stain = os.path.join(self.input, 'solid-ocean.out')
        self.stress_txt_output = os.path.join(self.stress_output, study_region + '_stress_Vol_N_S.txt')

        # Figure files
        self.select_data_fig = os.path.join(self.catalog_output, study_region + '_select_catalog.png')
        self.decluster_data_fig_lat = os.path.join(self.catalog_output, study_region + '_decluster_catalog_lat.png')
        self.decluster_data_fig_cum = os.path.join(self.catalog_output, study_region + '_decluster_catalog_cum.png')
        self.decluster_data_fig_nna = os.path.join(self.catalog_output, study_region + '_decluster_catalog_nna.png')

        self.all_region_vol_pdf = os.path.join(self.figure_output, study_region + '_all_region_Vol.pdf')
        self.all_region_S_pdf = os.path.join(self.figure_output, study_region + '_all_region_S.pdf')
        self.all_region_N_pdf = os.path.join(self.figure_output, study_region + '_all_region_N.pdf')
        self.all_region_CFS_pdf = os.path.join(self.figure_output, study_region + '_all_region_CFS.pdf')
        self.all_region_tidal_sensitivity_pdf = os.path.join(self.figure_output, study_region + '_all_region_tidal_sensitivity.pdf')

        self.all_region_tidal_sensitivity_vol_pdf = os.path.join(self.figure_output, study_region + '_all_region_tidal_sensitivity_vol.pdf')
        self.all_region_tidal_sensitivity_n_pdf = os.path.join(self.figure_output, study_region + '_all_region_tidal_sensitivity_n.pdf')
        self.all_region_tidal_sensitivity_s_pdf = os.path.join(self.figure_output, study_region + '_all_region_tidal_sensitivity_s.pdf')
        self.all_region_tidal_sensitivity_cfs_pdf = os.path.join(self.figure_output, study_region + '_all_region_tidal_sensitivity_cfs.pdf')

        self.time_tidal_P_value_pdf = os.path.join(self.figure_output, study_region + '_time_tidal_P_value.pdf')
        self.time_tidal_sensitivity_pdf = os.path.join(self.figure_output, study_region + '_time_tidal_sensitivity.pdf')
        self.time_tidal_phase_shift_pdf = os.path.join(self.figure_output, study_region + '_time_tidal_phase_shift.pdf')
        self.time_tidal_amplitude_pdf = os.path.join(self.figure_output, study_region + '_time_tidal_amplitude.pdf')
        self.time_earthquake_number_pdf = os.path.join(self.figure_output, study_region + '_time_earthquake_number.pdf')

        self.all_region_stress_pdf = os.path.join(self.figure_output, study_region + '_time_tidal_stress.pdf')
        self.time_all_CFS_pdf = os.path.join(self.figure_output, study_region + '_time_all_CFS.pdf')
        self.earthquake_b_positive_pdf = os.path.join(self.figure_output, study_region + '_earthquake_b_positive.pdf')
