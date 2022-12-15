
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import ast


def create_exfs(config_dir, ecco_dir, print_level):

    L1_model_name = 'L1_mac_delta'

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils','init_file_creation'))

    start_year = 1992
    start_month = 1

    final_year = 1992
    final_month = 1

    # # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    # import create_L1_exf_field_ref as ebcr
    # ebcr.create_L1_exf_ref_file(config_dir, L1_model_name, print_level)

    proc_ids = np.arange(10).tolist()

    import create_L1_monthly_exfs as cef
    for proc_id in proc_ids:
        cef.create_L1_exf_fields(config_dir, ecco_dir, L1_model_name, proc_id,
                                 start_year, final_year, start_month, final_month, print_level)


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

    create_exfs(config_dir, ecco_dir, print_level=5)
   

