# TiMEpy

This package provides tools for analyzing Tidal Modulation of slow and fast Earthquakes.

<img src="https://raw.githubusercontent.com/Weifan-Lu/TiMEpy/main/Logo.png" alt="TiMEpy Logo" width="400"/>


## Quick Start

Follow these steps to get started quickly:

1. **Clone the repository and install dependencies:**

  ```
       git clone https://github.com/Weifan-Lu/TiMEpy.git
       pip install numpy matplotlib scipy
       cd TiMEpy/ex_ridgecrest
  ```
2. **Run the complete analysis:**
   
   All scripts are integrated
   ```
      python main_run.py
   ```
3. **Run each step separately:**
   ### Preprocessing
   - **Create output:**  
     ```
        python ex_pre0_create_output.py
     ```
   - **Select catalog:**  
     ```
        python ex_pre1_select_catalog.py
     ```
   - **Decluster using NNA:**  
     ```
        python ex_pre2_decluster_NNA.py
     ```
   - **Convert strain to stress:**  
     ```bash 
        python ex_pre3_strain_to_stress.py
     ```
   ### Analysis
   - **Calculate tidal phase:**
      ```
        python ex_ana1_calc_tidal_phase.py
      ```
   - **Entire region analysis:**
      ```
         python ex_ana2_entire_region.py
      ```
   - **Temporal variation analysis:**
      ```
          python ex_ana3_temp_variation.py
      ```
   - **b-value computation:**
      ``` 
         python ex_ana4_b_value.py
      ```

## Additional Tutorial
 ### Ridgecrest example
  
For a comprehensive, step-by-step tutorial on using TiMEpy, refer to the [Simple_tutorial.ipynb](https://github.com/Weifan-Lu/TiMEpy/blob/main/Simple_tutorial.ipynb)  file—a detailed guide developed using Jupyter Notebook

 ### How to get tidal strain
 
In addition to the earthquake catalog, tidal strain is also essential. Tidal strain is generated using **TidalStrain.2**, developed by Dr. Fuyuki Hirose. Once you have completed testing the examples and are ready to analyze your actual data, please use **TidalStrain.2** to obtain the tidal strain. Detailed analysis and practical experience are documented within the **TidalStrain.2** file.


## Contact

For further information or support, please contact [luweifan001@gmail.com]
Homepage: https://weifan-lu.github.io/
