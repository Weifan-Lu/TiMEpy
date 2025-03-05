# TiMEpy

This package provides tools for analyzing Tidal Modulation of slow and fast Earthquakes.

<img src="https://raw.githubusercontent.com/Weifan-Lu/TiMEpy/main/Logo.png" alt="TiMEpy Logo" width="400"/>


## Quick Start

Follow these steps to get started quickly:

1. **Clone the repository and install dependencies:**

   ```bash
       git clone https://github.com/Weifan-Lu/TiMEpy.git
       pip install numpy matplotlib scipy
       cd TiMEpy/ex_ridgecrest

2. **Run the complete analysis:**
   ```bash
    python main_run.py

3. **Run each step separately:**
   ### Preprocessing
   - **Create output:**  
     ```bash
        python ex_pre0_create_output.py
     
   - **Select catalog:**  
     ```bash
        python ex_pre1_select_catalog.py
     
   - **Decluster using NNA:**  
     ```bash 
        python ex_pre2_decluster_NNA.py
     
   - **Convert strain to stress:**  
     ```bash 
        python ex_pre3_strain_to_stress.py
   ### Analysis
   - **Calculate tidal phase:**
      ```bash
        python ex_ana1_calc_tidal_phase.py
      
    python ex_ana2_entire_region.py
    python ex_ana3_temp_variation.py
    python ex_ana4_b_value.py
