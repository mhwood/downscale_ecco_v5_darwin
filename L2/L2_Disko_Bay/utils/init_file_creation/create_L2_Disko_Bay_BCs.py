
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys
import ast


def create_BCs(config_dir, L2_model_name, parent_model_level, parent_model_name, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'L2', 'utils','init_file_creation'))

    start_year = 1992
    start_month = 6
    start_day = 1

    final_year = 1992
    final_month = 12
    final_day = 31

    # step 2: using the reference dict, organize downscaled BC into daily files
    if parent_model_level=='L1':

        ############################################################################################
        # Create the BC fields (3 steps)

        # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
        import create_L2_BC_field_ref as ebcr
        ebcr.create_L2_BC_ref_file(config_dir, L2_model_name, parent_model_level, parent_model_name, print_level)

        # proc_ids = np.arange(27).tolist()
        # # proc_ids = np.arange(12).tolist()
        #
        # import create_L2_daily_bcs_from_L1_ref_grid as cef
        # for proc_id in proc_ids:  # 7
        #     cef.create_bc_fields_via_interpolation(config_dir, L2_model_name, parent_model_name, proc_id,
        #                                            start_year, final_year, start_month,
        #                                            final_month, start_day, final_day, print_level)
        #
        # # step 3: combine all of the BC fields into a single file
        # import combine_and_rotate_L2_daily_bc_files as com
        # for proc_id in proc_ids:
        #     com.combine_and_rotate_L2_daily_bcs(config_dir, L2_model_name, proc_id,
        #                                         start_year, final_year, print_level)

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
   

