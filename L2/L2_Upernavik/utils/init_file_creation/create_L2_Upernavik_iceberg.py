
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def read_grid_geometry_from_nc(config_dir,model_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    hFacC = np.array(ds.variables['HFacC'][:, :, :])
    delR = np.array(ds.variables['drF'][:])
    ds.close()
    return(XC,YC,Depth,hFacC,delR)

def create_L2_iceberg_files(config_dir, L2_model_name, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    sys.path.insert(1, os.path.join(config_dir,'L2', 'utils','init_file_creation'))

    XC, YC, Depth, hFacC, delR = read_grid_geometry_from_nc(config_dir, L2_model_name)


    min_row = 10
    max_row = 80
    min_col = 360
    max_col = 440

    iceberg_mask = np.zeros(XC.shape)
    iceberg_mask[min_row:max_row, min_col:max_col] = 1
    iceberg_mask[hFacC[0,:,:]==0] = 0  # Mask out land points
    iceberg_mask[Depth<50] = 0  # Mask out shallow points

    rows = np.arange(375)
    cols = np.arange(450)
    Cols, Rows = np.meshgrid(cols, rows)
    plt.pcolormesh(Cols, Rows, iceberg_mask)
    plt.contour(Cols, Rows, Depth, levels=[1], colors='k')
    plt.show()













if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_L3_iceplume_files(config_dir)
   

