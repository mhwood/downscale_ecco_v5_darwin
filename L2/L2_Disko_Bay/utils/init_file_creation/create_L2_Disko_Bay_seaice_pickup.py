
import os
import argparse
import sys
import ast

def create_seaice_pickup_file(config_dir, L2_model_name, parent_model_level, parent_model_name, parent_model_pickup_iteration, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L2','utils','init_file_creation'))

    ####################################################################################################################
    # this code reads the L1 pickup file as a grid
    if parent_model_level == 'L1':
        import create_L2_seaice_pickup_from_L1 as c31
        c31.create_L2_seaice_pickup_file(config_dir, parent_model_name, parent_model_pickup_iteration, L2_model_name, print_level)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_seaice_pickup_file(config_dir)
   

