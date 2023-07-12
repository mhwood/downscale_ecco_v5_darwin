
import os

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys
import netCDF4 as nc4

def create_bathy_file(config_dir, print_level):

    model_name = 'L2_Santa_Barbara'

    if print_level >= 1:
        print('Creating the bathymetry file for the '+model_name+' model')

    if print_level >= 1:
        print('    - Reading from the file subsetted from the GEBCO 15 arc sec grid')

    ds = nc4.Dataset(os.path.join(config_dir,'L2',model_name,'input',model_name+'_bathymetry_GEBCO.nc'))
    elev = ds.variables['elevation'][:, :]
    ds.close()

    elev = np.array(elev)

    output_file = os.path.join(config_dir, 'L2', model_name, 'input', model_name + '_bathymetry.bin')
    elev.ravel('C').astype('>f4').tofile(output_file)







if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_mitgrid(config_dir)
   

