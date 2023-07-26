
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys
import ast


def create_BCs(config_dir, ecco_dir, print_level):

    L1_model_name = 'L1_East_Pacific'

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils','init_file_creation'))

    start_year = 2000
    start_month = 1

    final_year = 2005
    final_month = 12

    proc_ids = np.arange(105).tolist()  # (31+4)*3 = 105
    boundaries = ['north', 'west','south']

    # # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    # import create_L1_BC_field_ref as ebcr
    # ebcr.create_L1_BC_ref_file(config_dir, L1_model_name, print_level)

    import create_L1_monthly_BCs as cef
    for proc_id in proc_ids:
        cef.create_bc_fields_via_interpolation(config_dir, ecco_dir, L1_model_name, boundaries, proc_id,
                                               start_year, final_year, start_month, final_month, print_level,
                                               write_to_unbalanced=False)

    # step 3: combine all of the BC fields into a single file
    import combine_and_rotate_L1_monthly_bc_files as com
    var_names = ['THETA', 'SALT', 'UVEL', 'VVEL']
    for i in range(1, 32):
        var_names.append('PTRACE' + '{:02d}'.format(i))
    com.combine_and_rotate_L1_monthly_bcs(config_dir, L1_model_name, boundaries, var_names,
                                          start_year, final_year, print_level,
                                          read_from_unbalanced=False)

    import balance_L1_vector_BCs as bbcs
    bbcs.balance_bc_fields(config_dir, L1_model_name, boundaries, start_year, final_year, print_level)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_dir", action="store",
                        help="The directory where ECCO files are stored.", dest="ecco_dir",
                        type=str, required=True)


    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_dir = args.ecco_dir

    create_BCs(config_dir, ecco_dir, print_level=4)
   

