
import os
import argparse
import sys
import ast
import numpy as np

def create_exf_files(config_dir, L2_model_name, parent_model_level, parent_model_name, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L2','utils','init_file_creation'))

    start_year = 1993
    final_year = 1993

    ###################################################################################
    # The exf fields are created in 3 steps

    if parent_model_level == 'L1':

        ############################################################################################
        # no need to make a list of files because we're not using the dv output from L1 (only L0)

        ############################################################################################
        # interpolate from annual L1 input files as a grid
        proc_ids = np.arange(10)
        proc_ids = [0,1,2,3,4,5,6,8,9]
        # proc_ids = [9]

        import create_L2_daily_exf_files_from_L1 as cef
        for proc_id in proc_ids:
            cef.create_exf_fields_via_interpolation(config_dir, L2_model_name, parent_model_name, proc_id,
                                                    start_year, final_year, print_level)

        # combine all of the exf fields into a single file
        import combine_and_rotate_L2_daily_exf_files as com
        for proc_id in proc_ids:
            com.combine_and_rotate_L2_daily_exf_files(config_dir, L2_model_name, proc_id, start_year, final_year, print_level)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

