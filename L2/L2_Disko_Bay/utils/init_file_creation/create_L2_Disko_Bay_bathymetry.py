
import os
import argparse
import sys
import numpy as np

def create_bathy_file(config_dir, print_level):

    sys.path.insert(1, os.path.join(config_dir,'utils','init_file_creation'))
    import create_bathymetry as cb

    model_name = 'L3_Disko_Bay'
    level_name = 'L3'

    # these will be used to find the part of the domain which is connected
    # all non-connected parts will be masked
    central_wet_row = 180
    central_wet_col = 300

    hFacMinDr = 1.0
    hFacMin = 0.3

    delR = np.array([1.00, 1.14, 1.30, 1.49, 1.70])

    cb.create_bathymetry_file(config_dir, level_name, model_name,
                              central_wet_row, central_wet_col, hFacMinDr, hFacMin, delR,
                              print_level)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_bathy_file(config_dir)
   

