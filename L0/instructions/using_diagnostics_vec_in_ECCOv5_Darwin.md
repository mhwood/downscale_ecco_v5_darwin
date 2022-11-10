# Using the diagnostics_vec package for in ECCOv5 Darwin

The steps required to run ECCOv5 Darwin with diagnostics_vec are as follows:
1. Build ECCOv5 Darwin with diagnostics_vec
2. Generate 


## Generate the diagnostic_vec masks
To generate masks for the boundaries and surface variables, use the `create_diagnostic_vec_masks.py` function in the L0/utils directory. This script is run from the command line:
```
python3 create_diagnostic_vec_masks.py -d [config_dir] -e [ecco_dir] -m [list of model names]
```
The three required inputs to this script are:
 - config_dir: the directory where the L0 and L1 directories are stored (the head of this repo if you clone it)
 - ecco_dir: a path to a directory where ecco tile files are stored (must contain tile001.mitgrid through tile005.mitgrid)
 - model names: the list of model names to create masks for (e.g. L1_W_Greenland, L1_GOM, etc)

If desired, plots of these masks can additionally be 
