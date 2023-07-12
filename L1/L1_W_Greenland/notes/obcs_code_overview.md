# Overview of Boundary Condition Creation

To generate the BCs, first, a dictionary is created to reference the files generated in the L0 run with their destination boundary condition files for the L1 model:
```
python3 init_file_creation/create_daily_obcs_ref.py -d /path/to/config/dir
```

Next, monthly boundary condition files are created using the dictionary as follows:
```
python3 init_file_creation/create_daily_bcs_from_ref.py -d /path/to/config/dir 
```

Once the monthly boundary conditions have been generated, they can be combined with the following script:
```
python3 init_file_creation/combine_and_rotate_monthly_obcs_files.py -d /path/to/config/dir 
```
Finally, we must make small adjustments to the boundary conditions on the boundaries to minimize mass flux into the domain: when interpolating the velocities from the L0 to the denser grid of the L1, there is no guarantee the fluxes across the boundary will be balance. The following code calculates the flux imbalance into the domain and removes the flux imbalance from the velocity fields:
```
python3 init_file_creation/balance_boundary_fluxes.py -d /path/to/config/dir 
```
