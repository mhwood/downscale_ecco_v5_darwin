
import os

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys
import netCDF4 as nc4

def create_exch2_file(config_dir, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import create_exch2_file as cex2

    model_name = 'L2_Santa_Barbara'
    level_name = 'L2'
    sNx = 60
    sNy = 30

    if print_level >= 1:
        print('Creating the exch2 file and blanklist for the ' + model_name + ' model')

    cex2.create_exch2_file(config_dir, level_name, model_name, sNx, sNy)








if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_mitgrid(config_dir)
   

