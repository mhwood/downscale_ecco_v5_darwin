
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys
import netCDF4 as nc4


def create_mitgrid(config_dir, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import create_mitgrid_from_XY as cm

    model_name = 'L2_Monterey_Bay'

    # pass this a grid on the subdomain in EPSG 3413 coordinates
    if print_level >= 1:
        print('Creating the mitgrid file for the '+model_name+' model')

        # pass this a grid on the subdomain in EPSG 3413 coordinates
        if print_level >= 1:
            print('Creating the mitgrid file for the ' + model_name + ' model')

        if print_level >= 1:
            print('    - Generating the grid in polar coordinates')

        ds = nc4.Dataset(os.path.join(config_dir, 'L2', model_name, 'input', model_name + '_bathymetry_GEBCO.nc'))
        x = ds.variables['lon'][:]
        y = ds.variables['lat'][:]
        ds.close()

        XC, YC = np.meshgrid(x, y)

        x_resolution = np.mean(np.diff(x))
        y_resolution = np.mean(np.diff(y))

        XG, YG = np.meshgrid(x - x_resolution / 2, y - y_resolution / 2)

        XG = np.hstack([XG, XG[:, -1:] + x_resolution])
        XG = np.vstack([XG, XG[-1:, :]])

        YG = np.hstack([YG, YG[:, -1:]])
        YG = np.vstack([YG, YG[-1:, :] + y_resolution])

        if print_level >= 1:
            print('    - The C grid has ' + str(np.shape(XC)[0]) + ' rows and ' + str(np.shape(XC)[1]) + ' cols')
            print(
                '        - x: ' + str(np.min(x)) + ' to ' + str(np.max(x)) + ' (resolution:' + str(x_resolution) + ')')
            print(
                '        - y: ' + str(np.min(y)) + ' to ' + str(np.max(y)) + ' (resolution:' + str(y_resolution) + ')')
            print('    - The G grid has ' + str(np.shape(XG)[0]) + ' rows and ' + str(np.shape(XG)[1]) + ' cols')

        cm.create_mitgrid_file(config_dir, 'L2', model_name, XC, YC, XG, YG, print_level)








if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_mitgrid(config_dir)
   

