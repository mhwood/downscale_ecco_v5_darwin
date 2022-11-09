# Building the ECCOv5 Darwin model with diagnostics_vec

These instructions are nearly identical to those described in the [ECCOv5 Darwin Github](https://github.com/MITgcm-contrib/ecco_darwin/tree/master/v05/llc270) - the only major change is that the [diagnostics_vec](https://github.com/mhwood/diagnostics_vec) package is added to output pertinent variables for the boundary of the downscaled domain. 

Here, the steps are outlined to run this model on the Pleiades computing cluster. If the model is to be run on different machines, small changes might be necessary. 

## Gathering the ECCO Darwin Code
To start, follow the first few lines of the step 1 instructions in the ECCOv5 Darwin Github to obtain the model files:
```
==============
# 1. Get code

git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout 24885b71
cd ..
```


## Adding `diagnostics_vec` to MITgcm
There are 2 main steps to add `diagnostics_vec` to the model

### Step 1: Add the `diagnostics_vec` package files to MITgcm
The `diagnostics_vec` package files can be easily added to MITgcm using a Python utility function in the `diagnostics_vec` directory:
```
git clone https://github.com/mhwood/diagnostics_vec.git
cd diagnostics_vec/utils/
python3 copy_pkg_files_to_MITgcm.py -m ../../darwin3
cd ../../darwin3
```
A full summary of the changes made to MITgcm src/model files by the above code is shown below.

### Step 2: Edit compile time files for diagnostics_vec
There are 3 compile time files to edit for diagnostics_vec

1. `packages.conf`: This file is inside the `../ecco_darwin/v05/llc270/code_darwin/` directory. For this file, simply add a line for `diagnostics_vec`.

2. `DIAGNOSTICS_VEC_SIZE.h`: This is a new file, provided in this repository, which we will copy from the default `pkg` directory to the `../ecco_darwin/v05/llc270/code` directory. Next, edit the file so that reflects the number of `diagnostics_vec` masks which will be used in the model run, e.g.:
```
cp pkg/diagnostics_vec/DIAGNOSTICS_VEC_SIZE.h ../ecco_darwin/v05/llc270/code
vim ../ecco_darwin/v05/llc270/code/DIAGNOSTICS_VEC_SIZE.h
# File contents after modifications for this example:
C------------------------------------------------------------------------------|
C                           DIAGNOSTICS_VEC_SIZE.h
C------------------------------------------------------------------------------|

      INTEGER, PARAMETER :: VEC_points = 900
      INTEGER, PARAMETER :: nVEC_mask = 3
      INTEGER, PARAMETER :: nSURF_mask = 1
```

3. `CPP_OPTIONS.h`: Edit this file if external forcing fields are desired for the output. This file is inside the `../ecco_darwin/v05/llc270/code` directory. For this file, we will add three lines so that `diagnostics_vec` can access external forcing variables (not defined when the `ecco` package is used for some reason I can't figure out...):
```
18 C-- Forcing code options:
19 
20 #define ALLOW_ATM_TEMP                                 # added line
21 #define ALLOW_DOWNWARD_RADIATION                       # added line
22 #define ALLOW_RUNOFF                                   # added line
```

## Finish the model build

After the package is added and code modification files are edited, the model can be rebuilt using the same commands in ECCOv5 Darin [readme](https://github.com/MITgcm-contrib/ecco_darwin/blob/master/v05/llc270/readme.txt) which are copied here for convenience:
```
 ==============
 # 2. Build executable
 
 # should be inside the darwin3 dir
 mkdir build run
 cd build

 module purge
 module load comp-intel mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
 ../tools/genmake2 -of ../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas \
   -mo '../../ecco_darwin/v05/llc270/code_darwin ../../ecco_darwin/v05/llc270/code' -mpi
 make depend
 make -j 16
 
 ```
 
## Run the model
Now that the model is built, it is nearly ready to run. To run this model, follow the steps provided in the running_ECCOv5_darwin_with_dv.md instructions.


---


### Summary of MITgcm Changes for `diagnostics_vec` (

1. `PARAMS.h`: This file is inside the `model/inc` directory. Here, two lines are added for the `diagnostics_vec` package. First `useDiagnostics_vec` is defined as a LOGICAL and then put it in the namelist:
```
1057      LOGICAL useDiagnostics
1058      LOGICAL useDiagnostics_vec                                  # new line added
1059      LOGICAL useREGRID
```
and
```
1076     &        useDiagnostics, useREGRID, useLayers, useMNC,
1077     &        useDiagnostics_vec,                                  # new line added
1078     &        useRunClock, useEMBED_FILES,
```

2. `packages_boot.F`: This file is inside the `model/src` directory. For this file, the `diagnostics_vec` package is added in the same location where the `diagnostics` package is included. This occurs in three places:
```
91     &          useDiagnostics,
92     &          useDiagnostics_vec,                  # added line
93     &          useREGRID,
```
and
```
157      useDiagnostics  =.FALSE.
158      useDiagnostics_vec  =.FALSE.                         # added line
159      useREGRID       =.FALSE.
```
and
```
393 #ifdef ALLOW_DIAGNOSTICS
394       CALL PACKAGES_PRINT_MSG( useDiagnostics,'Diagnostics', ' ' )
395 #endif
396 #ifdef ALLOW_DIAGNOSTICS_VEC                                                   # added line
397       CALL PACKAGES_PRINT_MSG( useDiagnostics_vec,                             # added line
398      &                         'Diagnostics_vec', ' ' )                        # added line
399 #endif                                                                         # added line
400 #ifdef ALLOW_REGRID
401       CALL PACKAGES_PRINT_MSG( useREGRID,     'REGRID',      ' ' )
402 #endif
```

3. `packages_init_fixed.F`: This file is inside the `model/src` directory. For this file, the `diagnostics_vec` package is added near the end of the file before the `ALLOW_DIAGNOSTICS` block:
```
662 #endif /* ALLOW_CTRL */
663 
664 #ifdef ALLOW_DIAGNOSTICS_VEC                                                   # added line
665       IF ( useDiagnostics_vec ) THEN                                           # added line
666 # ifdef ALLOW_DEBUG                                                            # added line
667         IF (debugMode)                                                         # added line
668      & CALL DEBUG_CALL('DIAGNOSTICS_VEC_INIT_FIXED',myThid)                    # added line
669 # endif                                                                        # added line
670         CALL DIAGNOSTICS_VEC_INIT_FIXED( myThid )                              # added line
671       ENDIF                                                                    # added line
672 #endif                                                                         # added line
673 
674 #ifdef ALLOW_DIAGNOSTICS
```

4. `packages_readparms.F`: This file is inside the `model/src` directory. For this file, the `diagnostics_vec` package is added fter the block for the `diagnostics` package:
```
370 #endif /* ALLOW_DIAGNOSTICS */
371 
372 #ifdef ALLOW_DIAGNOSTICS_VEC
373       CALL DIAGNOSTICS_VEC_READPARMS( myThid )
374 #endif /* ALLOW_DIAGNOSTICS_VEC */
```

5. `packages_check.F`: This file is inside the `model/src` directory. For this file, the `diagnostics_vec` package is added after the block for the `diagnostics` package:
```
434 #ifdef ALLOW_DIAGNOSTICS_VEC
435       IF (useDiagnostics_vec) CALL DIAGNOSTICS_VEC_CHECK( myThid )
436 #else
437       IF (useDiagnostics_vec)
438      &   CALL PACKAGES_ERROR_MSG( 'Diagnostics_vec', ' ', myThid )
439 #endif
```

6. `do_the_model_io.F`: This file is inside the `model/src` directory. For this file, the `diagnostics_vec` package is added before the block for the `diagnostics` package:
```
268 #ifdef ALLOW_DIAGNOSTICS_VEC
269       IF ( useDiagnostics_vec )
270      &     CALL DIAGNOSTICS_VEC_OUTPUT( myTime, myIter, myThid )
271 #endif
```

