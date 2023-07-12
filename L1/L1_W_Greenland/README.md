# The L1_W_Greenland model

This configuration is constructed to simulate ocean circulation and biogeochemistry in Davis Strait and Baffin Bay, West of Greenland.

The following steps will set up the grid, bathymetry, initial conditions, boundary conditions, and external forcing conditions for the model run. Note that the [create_all_L1_W_Greenland_files.py](https://github.com/mhwood/downscale_ecco_v5_darwin/blob/main/L1/L1_W_Greenland/utils/init_file_creation/create_all_L1_W_Greenland_files.py) script provides a pipeline to run all of the pertinent scripts for the files for this model.

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
The initial condtions, boundary conditions, and external forcing conditions all must be interpolated from the "parent" model to the new downscaled model. To easily generate the grid as a complete netcdf file, my strategy is to set up the model and run for one time step, outputting the grid using the mnc package. For this model run, 

### Step 3: Create the initial conditions
The initial conditions for the subdomain are stored in a pickup file for convenience. To generate this pickup file, an equivalent pickup file from the L0_540 run is interpolated into the new grid. For our investigation, we used a pickup file 1 year (175,680 time steps) into the L0_540 simulation, located in the L0_540/run directory: pickup.0000175680.data. Note also that the code will use the tile004.mitgrid and tile005.mitgrid files, which are expected to be located in the L0_540/run directory. To generate the new pickup file, run the following code:
```
python3 init_file_creation/create_pickup_file.py -d /path/to/config/dir -i 175680
```
This script will generate the files ```pickup.0000351360.data``` and ```pickup.0000351360.meta```. Note that the iteration counter has been doubled because the time step of the L1 model is half that of the L0 model. If desired, a plot showing all of the initial condition fields can be made as follows:
```
python3 plot_creation/plot_pickup_file_components.py -d /path/to/config/dir -i 351360
```

### Step 4: Create the external forcing fields
The external forcing fields for this configuration are downscaled from output from the L0 model run. After running the L0 model, th `L0_540/run` directory will have a number of files containing approximately monthly external forcing fields, which were generated by the ```diagnostics_vec``` package. First, we create a dictionary that references which files will be used to contruct the monthly external forcing conditions for the L1 model. To generate this dictionary, run the following script:
```
python3 init_file_creation/create_monthly_exf_field_ref.py -d /path/to/config/dir
```
Note that the dictionary (called exf_dest_ref.txt) is created in the ```L0/run/dv``` directory. Next, the dictionary is used to generate monthly external forcing fields by specifying the time bounds of interest:
```
python3 init_file_creation/create_monthly_exf_fields_from_ref.py -d /path/to/config/dir -i [field number (0-10)] -S [start year] -s [start month] -F [end year] -f [end month]
```

Note that the script is created to generate one field at a time over a specified duration. The script was written in this way so that it could be run in parallel on a high end computing cluster. Further, a date range is explicitly specified. For our investigation, we ran the L1 model from a start year/month of 1993/1 and an end year/month of 1996/12, i.e to create the atmospheric temperature conditions, use:
```
python3 init_file_creation/create_monthly_exf_fields_from_ref.py -d /path/to/config/dir -i 2 -S [start year] -s [start month] -F [end year] -f [end month]
```
Once the monthly external fields have been generated, they can be combined with the following script:
```
python3 init_file_creation/combine_monthly_exf_files.py -d /path/to/config/dir -i [field number (0-10)] -S [start year] -s [start month] -F [end year] -f [end month]
```
If desired, a plot showing an external forcing fields at a particular timestep can be generated using the following function:
```
python3 plot_creation/plot_exf_field.py -d /path/to/config/dir -i [field number (0-10)] -t [time_step]
```


### Step 5: Create the boundary conditions
The creation of the boundary conditions for this configuration follow a similar set of steps used to generate the external forcing conditions. First, a dictionary is created to reference link the files generated in the L0_540 run with their destination boundary condition files for the L1_1080 model:
```
python3 init_file_creation/create_daily_obcs_ref.py -d /path/to/config/dir
```
Next, daily boundary condition files are created using the dictionary as follows:
```
python3 init_file_creation/create_daily_bcs_from_ref.py -d /path/to/config/dir -i [field number (0-15)] -S [start year] -s [start month] -sd [start day] -F [end year] -f [end month] -fd [end day]
```
The list of field numbers correspond to the following variables/masks:

0. THETA on the northern boundary
1. SALT on the southern boundary
2. UVEL on the western boundary
3. VVEL on the northern boundary
4. ETAN on the southern boundary   (not used)
5. THETA on the western boundary
6. SALT on the northern boundary
7. UVEL on the southern boundary
8. VVEL on the western boundary
9. ETAN on the northern boundary   (not used)
10. THETA on the southern boundary
11. SALT on the western boundary
12. UVEL on the northern boundary
13. VVEL on the southern boundary
14. ETAN on the western boundary   (not used)

For example, to general the meridional velocity across the northern boundary for the entirety of the 4-year model time domain, run the following line
```
python3 init_file_creation/create_daily_bcs_from_ref.py -d /path/to/config/dir -i 3 -S 1993 -s 1 -sd 1 -F 1996 -f 12 -fd 31
```
Once the daily boundary conditions have been generated, they can be combined with the following script:
```
python3 init_file_creation/combine_daily_obcs_files.py -d /path/to/config/dir -i [field number (0-15)] -S [start year] -s [start month] -s [start day] -F [end year] -f [end month] -fd [end day]
```
Finally, we must make small adjustments to the boundary conditions on the west, north and south boundaries. When inteprolating the velocities from the L0_540 to the denser grid of the L1_1080, there is no guarantee the fluxes across the boundary will be balance. In addition, we cannot use the online MITgcm code to balance boundary fluxes, because our model has barotropic tides. Thus, we use the following code, which calculates the flux imbalance into the domain, smooths it weekly, and then removes the smooth flux imbalance from the velocity fields:
```
python3 init_file_creation/balance_boundary_fluxes.py -d ../../
```

If desired, a plot showing the boundary conditions at a particular timestep can be generated using the following function:
```
python3 plot_creation/plot_boundary_conditions.py -t [time_step]
```
Note that the boundary conditions are averages in 15 minute internals - the temporal interpolation, similar to the external forcing fields, is accomodated by the obcs package.

### Step 6: Create the diagnostics_vec files
In order to configuration the subsequent model level (L2), the boundary conditions are generated using the diagnostics_vec package:
```
python3 init_file_creation/generate_diagnostics_vec_masks.py
```