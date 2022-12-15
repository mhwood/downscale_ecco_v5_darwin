
import os
import netCDF4 as nc4
import numpy as np
import argparse
import sys

def plot_darwin_pickup(config_dir, print_level):

    L1_model_name = 'L1_mac_delta'

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils','plot_creation','init_files'))

    pickup_iteration = 292

    # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    import plot_L1_darwin_pickup_field as pf
    pf.create_darwin_pickup_plot(config_dir, L1_model_name, pickup_iteration)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_darwin_pickup(config_dir, print_level=4)
   

