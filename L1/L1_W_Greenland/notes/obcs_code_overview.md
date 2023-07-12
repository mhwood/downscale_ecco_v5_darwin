# Overview of Boundary Condition Creation

First, a dictionary is created to reference link the files generated in the L0_540 run with their destination boundary condition files for the L1_1080 model:
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
