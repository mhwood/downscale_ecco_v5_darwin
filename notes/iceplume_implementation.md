# iceplume implementation

This note is to describe manual changes to the src/inc files of the darwin3 fork of MITgcm to implement the iceplume package, which is currently activated in the L2 model. I have added the iceplume package manually here (rather than using An Nguyen's pull request fork) because a) I have already altered thee boot sequence files with `diagnostics_vec`, which is not yet checked in to MITgcm, b) I would like to understand how the package is implemented, and c) I made some modifications to the code after running into some complime-time errors.

The first several steps are just to modify the standard boot sequence. The notes at the end are specific changes to other scripts.

### Step 1: Add package files to the MITgcm/pkg directory

### Step 2: Add key words for `iceplume` in `PARAMS.h`
Define the key word:
```
       LOGICAL useICEPLUME
```
and add it to the `/PARM_PACKAGES/` namelist:
```
     &        useIceplume,
```

### Step 3: Add  `iceplume` into `packages_boot.F`
First add the keyword to the namelist,
```
     &        useIceplume,
```
initialize the keyword as false,
```
      useICEPLUME     =.FALSE.
```
and then add a block to call the package:
```
#ifdef ALLOW_ICEPLUME
      CALL PACKAGES_PRINT_MSG( useICEPLUME,  'ICEPLUME',   ' ' )
#endif
```

### Step 4: Add a block for ```iceplume``` into `packages_check.F`
Add this block:
```
#ifdef ALLOW_ICEPLUME
      IF (useICEPLUME) CALL ICEPLUME_CHECK( myThid )
#else
      IF (useICEPLUME) CALL PACKAGES_ERROR_MSG('ICEPLUME',' ',myThid)
#endif
```

### Step 5: Add a block for ```iceplume``` into `packages_init_fixed.F`
Add this block:
```
#ifdef ALLOW_ICEPLUME
      IF (useICEPLUME) THEN
# ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL('ICEPLUME_INIT_FIXED',myThid)
# endif
        CALL ICEPLUME_INIT_FIXED(myThid)
      ENDIF
#endif
```

### Step 6: Add a block for ```iceplume``` into `packages_init_variables.F`
Add this block:
```
#ifdef ALLOW_ICEPLUME
      CALL ICEPLUME_INIT_VARIA( myThid )
#endif
```

### Step 7: Add a block for ```iceplume``` into `packages_readparms.F`
Add this block:
```
#ifdef ALLOW_ICEPLUME
C--   if useICEPLUME=T, set mypackage parameters; otherwise just return
      CALL ICEPLUME_READPARMS( myThid )
#endif
```

### Step 8: Add the following `src` list of files from An Nguyen's modifications:
  - apply_forcing.F
  - do_oceanic_phys.F
  - external_fields_load.F
  - external_forcing.F
  - ini_parms.F

 The above files work well with one exception: on my machine the `temp_addMass3Dplume` and `salt_addMass3Dplume` cause compile time errors because they are too long. I have modified these key words to be `t_addMass3Dplume` and `s_addMass3Dplume`. These changes take place in the following files:
  - apply_forcing.F
  - ICEPLUME.h
  - iceplume_calc.F
  - iceplume_init_varia.F

### Step 9: Add the following `pkg` list of files from An Nguyen's modifications:
  - diagnostics/diagnostics_fill_state.F
  - diagnostics/diagnostics_main_init.F
  - exf/exf_getffields.F

### Final notes:
There are a few snags I ran into:
1. Ensure that `#ALLOW_ADDFLUID` is defined in `CPP_OPTIONS.h`
2. Ensure that `selectAddFluid = 1,` in the `data` file.
3. Make sure `exf` is turned on - it obviously is here but this is not always the case (e.g. when testing toy models).
4. There were some issues with the overlap in iceplume_calc.F. I edited some of the loops to includes the overlaps.
5. In `iceplume_readparms.F`, an `#ALLOW_EXCH2` header block needed to be moved above a `#USE_EXF_INTERPOLATION` block.
6. The pickup files need an `AddMass` field. If the model is initialized in winter, an AddMass field can be added which is all 0's.

