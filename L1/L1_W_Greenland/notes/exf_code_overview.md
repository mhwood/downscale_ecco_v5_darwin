# Overview of External Forcing Pipline

These scripts are called in succession in [create_L1_W_Greenland_exf.py](https://github.com/mhwood/downscale_ecco_v5_darwin/blob/main/L1/L1_W_Greenland/utils/init_file_creation/create_L1_W_Greenland_exf.py).

First, we create a dictionary that references which files will be used to contruct the monthly external forcing conditions for the L1 model. To generate this dictionary, run the following script:
```
python3 init_file_creation/create_L1_exf_ref_file.py -d /path/to/config/dir
```
Note that the dictionary (called exf_dest_ref.txt) is created in the ```L0/run/dv``` directory. Next, the dictionary is used to generate monthly external forcing fields by specifying the time bounds of interest:
```
python3 init_file_creation/create_L1_monthly_exfs.py -d /path/to/config/dir
```
Be sure to specify the years of interest.

The above script will create the exf conditions in monthly files, but MITgcm expects them in annual files for our configuration. To combine them into annual files, run:
```
python3 init_file_creation/combine_and_rotate_L1_monthly_exf_files.py -d /path/to/config/dir 
```
