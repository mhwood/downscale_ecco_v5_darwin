
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys
import ast


def create_ptracer_pickup(config_dir, ecco_dir, print_level):

    L1_model_name = 'L1_W_Greenland'

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils','init_file_creation'))

    parent_model_pickup_iteration = 73

    llc = 270
    ordered_ecco_tiles = [[101, 110], [100, 109], [62, 61]]
    ordered_ecco_tile_rotations = [[1, 1], [1, 1], [2, 2]]  # rotations are counter-clockwise

    import create_L1_ptracer_pickup as cp
    cp.create_L1_ptracer_pickup_file(config_dir, L1_model_name,
                               ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                               parent_model_pickup_iteration, print_level)


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

    create_ptracer_pickup(config_dir, ecco_dir, print_level=4)
   

