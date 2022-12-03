
import os
import netCDF4 as nc4
import numpy as np
import argparse
import sys

def plot_BCs(config_dir, print_level):

    L1_model_name = 'L1_East_Pacific'

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils','plot_creation','init_files'))

    year = 1992
    month = 1
    day = 17

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',L1_model_name+'_grid.nc'))
    hFac = ds.variables['HFacC'][:,:,:]
    ds.close()
    Nr = np.shape(hFac)[0]
    n_rows = np.shape(hFac)[1]
    n_cols = np.shape(hFac)[2]

    boundaries = ['north','west','south']

    # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    import plot_L1_BC_fields as pbc
    pbc.plot_L1_BCs(config_dir, L1_model_name, boundaries, Nr, n_rows, n_cols, year, month, day)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_BCs(config_dir, print_level=4)
   

