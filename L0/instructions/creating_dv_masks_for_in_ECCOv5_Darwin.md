# Creating diagnostic_vec masks for use in the ECCOv5 Darwin model


## Prepare subdomain grids
To generate the masks for model boundary conditions, the model grid is required. It is convenient to generate similar grid files for each subdomain model so that the mask script can easily access the subdomain geometries. To generate these files, I have found it is easier to run each model for one timestep, leveraging the nice `mnc` package that groups model grid parameters into a single file. Instructions for running each model to generate the grid are provided in each L1 model subdirectory.

Once the grids for each of these models is generate, a copy is made in a `downscale_darwin/nc_grids` directory.


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
