
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys
import ast


def create_BCs(config_dir, L1_model_name, L2_model_name, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'L2', 'utils','init_file_creation'))

    start_year = 1992
    start_month = 1

    final_year = 1993
    final_month = 12

    ############################################################################################
    # Create the BC fields (3 steps)

    # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    # import create_L2_BC_field_ref as ebcr
    # ebcr.create_L2_BC_ref_file(config_dir, L1_model_name, L2_model_name, print_level)

    proc_ids = np.arange(120).tolist()  # (31+9)*3 = 120
    boundaries = ['north', 'south', 'west']

    # import create_L2_monthly_BCs as cef
    # for proc_id in proc_ids:
    #     cef.create_bc_fields_via_interpolation(config_dir, L1_model_name, L2_model_name, boundaries, proc_id,
    #                                            start_year, final_year, start_month, final_month, print_level,
    #                                            write_to_unbalanced=False)

    # # step 3: combine all of the BC fields into a single file
    # import combine_and_rotate_L2_monthly_bc_files as com
    # var_names = ['THETA', 'SALT', 'UVEL', 'VVEL', 'UICE', 'VICE', 'HSNOW', 'HEFF', 'AREA']
    # for i in range(1, 32):
    #     var_names.append('PTRACE' + '{:02d}'.format(i))
    # com.combine_and_rotate_L2_monthly_bcs(config_dir, L2_model_name, boundaries, var_names,
    #                                       start_year, final_year, print_level,
    #                                       read_from_unbalanced=False)

    import balance_L2_vector_BCs as bbcs
    bbcs.balance_bc_fields(config_dir, L2_model_name, boundaries, start_year, final_year, print_level)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-i", "--proc_id", action="store",
                        help="The id of the process to run.", dest="proc_id",
                        type=int, required=True)

    parser.add_argument("-S", "--start_year", action="store",
                        help="The start year.", dest="start_year",
                        type=int, required=True)

    parser.add_argument("-s", "--start_month", action="store",
                        help="The start month.", dest="start_month",
                        type=int, required=True)

    parser.add_argument("-sd", "--start_day", action="store",
                        help="The start day.", dest="start_day",
                        type=int, required=True)

    parser.add_argument("-F", "--final_year", action="store",
                        help="The final year.", dest="final_year",
                        type=int, required=True)

    parser.add_argument("-f", "--final_month", action="store",
                        help="The final ymonth.", dest="final_month",
                        type=int, required=True)

    parser.add_argument("-fd", "--final_day", action="store",
                        help="The final day.", dest="final_day",
                        type=int, required=True)


    args = parser.parse_args()
    config_dir = args.config_dir
    proc_id = args.proc_id
    start_year = args.start_year
    final_year = args.final_year
    start_month = args.start_month
    final_month = args.final_month
    start_day = args.start_day
    final_day = args.final_day

    create_BCs(config_dir,proc_id,
               start_year, final_year, start_month,
               final_month, start_day, final_day)
   

