# Running the ECCOv5 Darwin model with diagnostics_vec

To run the ECCOv5 Darwin model with the diagnostics_vec package, begin by copying/linking all of the input binaries and runtime files as described on the ECCOv5 Darwin Github:
```
==============
# 3. Instructions for running simulation (1992-2020 period)

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp19/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_initial_conditions/* .
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_forcing/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2021/xx*42.data .
cp ../../ecco_darwin/v05/llc270/input/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
```

Next, we'll modify the run directory for the diagnostics_vec package in three steps:

### Step 1: Make masks and add them to a run/dv dir
The `diagnostics_vec` pkg ingests 2D masks which delineate where the model output is requested. The masks are stored in compact format, identical to the format of the bathymetry file. The masks are identically zero except where output is requested, in which case the mask is numbered sequentially (e.g. along a boundary). 

### Step 2: Add the data.diagnostics_vec file
After the masks are constructed, add the data.diagnostics_file to the run directory. This file lists the masks, the number of iterations to be stored in each output file, the averaging period for each mask, and the variables to output on each mask. Masks are split into two categories: lateral masks (e.g. for boundary conditions) which can output 2D or 3D variables, and surface masks which can access only variables at the surface (2D variables or the surface cells of 3D variables). An example diagnostics_vec file is shown here:
```
&DIAG_VEC_INPUT_VARS
#
 nml_startTime = 0,
 nml_endTime = 3153600000.,
#
# lateral BC's are averaged hourly and dumped monthly
#
 nml_vecFiles(1)  = 'dv/L1_W_north_mask.bin',
 nml_vecFiles(2)  = 'dv/L1_W_east_mask.bin',
 nml_vecFiles(3)  = 'dv/L1_W_south_mask.bin',
#
# lateral bc's have 720 iterations per file (24*30)
#
 nml_vec_iters_per_file(1:3) = 720, 720, 720,
#
# lateral bc's have an averaging period of 60 min (60*60)
#
 nml_vec_avg_periods(1:3) = 3600., 3600., 3600.,
#
 nml_fields2D(1,1)  = 'ETAN    ',
 nml_fields2D(1,2)  = 'ETAN    ',
 nml_fields2D(1,3)  = 'ETAN    ',
#
 nml_fields3D(1:2,1)  = 'PTRACE01','PTRACE02',
 nml_levels3D(1:2,1)  =   50, 50,
#
 nml_fields3D(1:2,2)  = 'PTRACE01','PTRACE02',
 nml_levels3D(1:2,2)  =   50, 50,
#
 nml_fields3D(1:2,3)  = 'PTRACE01','PTRACE02',
 nml_levels3D(1:2,3)  =   50, 50,
#
# surface bc's are averaged every 6 hour and dumped every month
#
 nml_sfFiles(1) = 'dv/L1_W_surface_mask.bin',
#
# surface bc's have 120 iterations per file (30*4)
#
 nml_sf_iters_per_file(1) = 120,
#
# surface bc's have an averaging period of 6 hour (6*60*60)
#
 nml_sf_avg_periods(1) = 21600.,
#
 nml_fieldsSF(1,1) = 'PCO2    ',
#
 nml_filePrec = 32,
 &
```

### Step 3: Add diagnostics_vec to data.pkg
As a final step, ensure that the diagnostics_vec package is activated in the data.pkg file

## Run the mode
After these three run time adjustments are made, the model is ready to run. Follow the usual job submission for Pleiades.
