# The L1_W_Greenland model

This configuration is constructed to simulate ocean circulation and biogeochemistry in Davis Strait and Baffin Bay, West of Greenland.

The following steps will set up the grid, bathymetry, initial conditions, boundary conditions, and external forcing conditions for the model run. Note that the [create_all_L1_W_Greenland_files.py](https://github.com/mhwood/downscale_ecco_v5_darwin/blob/main/L1/L1_W_Greenland/utils/init_file_creation/create_all_L1_W_Greenland_files.py) script provides a pipeline to run all of the pertinent scripts for the files for this model.

The scripts below will require two main file paths:
| Flag          | Expected Path                                                         |
|---------------|-----------------------------------------------------------------------|
| -e ecco_dir   | path to files from the LLC270 and LLC1080 models                      |
| -d config_dir | path to this repository (where the L0 and L1 directories are located) |

## Step 0: Prerequisite
Before the L1 configuration can be prepared, the L0 configuration must be run to output boundary and external forcing conditions. This example assumes the [diagnostics_vec](https://github.com/mhwood/diagnostics_vec) package has been used to generate daily boundary conditions.

## Step 1: Creating the mitgrid file
This configuration uses an mitgrid file to generate the model grid. To generate this file, use the following script:
file:
```
python3 init_file_creation/create_L1_W_Greenland_mitgrid.py -d /path/to/config/dir -e /path/to/ECCO/dir
```

## Step 2: Creating the bathymetry
This configuration uses bathymetry from the LLC1080 grid bathymetry. To generate this file, use the following script:
file:
```
python3 init_file_creation/create_L1_W_Greenland_bathymetry.py -d /path/to/config/dir -e /path/to/ECCO/dir
```

## Step 2.5 Intermediate Step: Make an nc file for the model grid
The initial condtions, boundary conditions, and external forcing conditions all must be interpolated from the "parent" model to the new downscaled model. To easily generate the grid as a complete netcdf file, my strategy is to set up the model and run for one time step, outputting the grid using the mnc package. For this model run, I use a paired-down version of the model using the [code_for_grid](https://github.com/mhwood/downscale_ecco_v5_darwin/tree/main/L1/L1_W_Greenland/code_for_grid) and [namelist_for_grid](https://github.com/mhwood/downscale_ecco_v5_darwin/tree/main/L1/L1_W_Greenland/namelist_for_grid) along with the input/L1_W_Greenland.mitgrid and input/L1_W_Greenland_bathymetry.bin files generated in the previous 2 steps. 

To run this model for one time step on 12 cpus (for example), split the mitgrid file into 180 by 180 tiles using:
```
python3 init_file_creation/split_L1_W_Greenland_mitgrid_file_to_tiles.py -d /path/to/config/dir -f 1 -r 180 -c 180
```
Then, build and run the model with the [build_for_grid](https://github.com/mhwood/downscale_ecco_v5_darwin/blob/main/L1/L1_W_Greenland/utils/build_for_grid.sh) and [run_for_grid](https://github.com/mhwood/downscale_ecco_v5_darwin/blob/main/L1/L1_W_Greenland/utils/run_for_grid.sh) scripts.


This model run will output pieces of the grid as nc files in the mnc* directories. To stitch these grids together into a netcdf file that spans the whole domain, use
```
python3 init_file_creation/stitch_L1_W_Greenland_nc_grid_files_for_ref.py -d /path/to/config/dir
```
This nc file will output in a general `nc_files` directory in the main configuration directory.


## Step 3.1: Create the initial conditions
The initial conditions for the subdomain are stored in a pickup file for convenience. To generate this pickup file, an equivalent pickup file from the ECCO_Darwin run is interpolated into the new grid. For our investigation, we used a pickup file a little more than 24 hours (73 time steps) into the ECCO_Darwin simulation, located in the L0/run directory: pickup.0000000073.data. To generate the new pickup file, run the following code:
```
python3 init_file_creation/create_pickup_file.py -d /path/to/config/dir -e /path/to/ECCO/dir
```
This script will generate the files ```pickup.0000000292.data``` and ```pickup.0000000292.meta```. Note that the iteration counter has been quadrupled because the time step of the L1 model is 4 times that of the L0 model. If desired, a plot showing all of the initial condition fields (at the surface) can be made as follows:
```
python3 plot_creation/init_files/plot_L1_W_Greenland_pickup_fields.py -d /path/to/config/dir
```

## Step 3.2: Create the seaice initial conditions
In a similar sense, the seaice pickup can be generated with:
```
python3 init_file_creation/create_seaice_pickup_file.py -d /path/to/config/dir -e /path/to/ECCO/dir
```

## Step 3.3: (Biogeochemistry Only) Create the biogeochemistry initial conditions
Similarly, the darwin and ptracers pickups can be generated with:
```
python3 init_file_creation/create_darwin_pickup_file.py -d /path/to/config/dir -e /path/to/ECCO/dir
python3 init_file_creation/create_ptracers_pickup_file.py -d /path/to/config/dir -e /path/to/ECCO/dir
```

## Step 4: Create the external forcing fields
The external forcing fields for this configuration are downscaled from output from the L0 model run. After running the L0 model, the `L0/run` directory will have a number of files containing approximately monthly external forcing fields, which were generated by the ```diagnostics_vec``` package. This output is interpolated to generate external forcing conditions for the downscaled model using the following code:
```
python3 init_file_creation/create_L1_W_Greenland_exf.py -d /path/to/config/dir -e /path/to/ECCO/dir
```
Be sure to edit this script to generate the correct conditions for your simulation. For physics-only runs, the external forcing fields are AQH, ATEMP, LWDOWN, PRECIP, RUNOFF, SWDOWN, UWIND, and VWIND; biogeochemistry adds IRONDUST and ATMOSCO2. 

The steps in this script are summarized in the [exf_code_overview]. (https://github.com/mhwood/downscale_ecco_v5_darwin/blob/main/L1/L1_W_Greenland/notes/exf_code_overview.md) file.

If desired, a plot showing an external forcing fields at a particular timestep can be generated using the following function:
```
python3 plot_creation/plot_exf_field.py -d /path/to/config/dir
```


## Step 5: Create the boundary conditions
The creation of the boundary conditions for this configuration follow a similar set of steps used to generate the external forcing conditions. 
```
python3 init_file_creation/create_L1_W_Greenland_BCs.py -d /path/to/config/dir -e /path/to/ECCO/dir
```
Be sure to edit this script to generate the correct conditions for your simulation. For physics-only runs, the boundary fields are THETA, SALT, UVEL, VVEL, AREA, HEFF, HSNOW, UICE, and VICE; biogeochemistry adds PTRACER01 through PTRACER31. 

The steps in this script are summarized in the [obcs_code_overview](https://github.com/mhwood/downscale_ecco_v5_darwin/blob/main/L1/L1_W_Greenland/notes/obcs_code_overview.md) file.

If desired, a plot showing an boundary conditions at a particular timestep can be generated using the following function:
```
python3 plot_creation/plot_BC_fields.py -d /path/to/config/dir
```

## Step 6: Create the diagnostics_vec files
In order to generate BCs for the subsequent model level (L2), the boundary conditions are generated using the diagnostics_vec package:
```
python3 init_file_creation/generate_diagnostics_vec_masks.py
```
These masks are used in the L1 run with the `diagnostics_vec` package.
