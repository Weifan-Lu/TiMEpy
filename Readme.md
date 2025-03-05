# TiMEpy

This package provides tools for analyzing Tidal Modulation of slow and fast Earthquakes.

![TiMEpy Logo]((https://github.com/Weifan-Lu/TiMEpy/blob/main/Logo.png)

    git clone https://github.com/Weifan-Lu/TiMEpy
    pip install numpy matplotlib scipy
    cd TiMEpy/ex_ridgecrest

Run all .py

    python main_run.py

Run each step separately

    python ex_pre0_create_output.py
    python ex_pre1_select_catalog.py
    python ex_pre2_decluster_NNA.py
    python ex_pre3_strain_to_stress.py
    python ex_ana1_calc_tidal_phase.py
    python ex_ana2_entire_region.py
    python ex_ana3_temp_variation.py
    python ex_ana4_b_value.py
