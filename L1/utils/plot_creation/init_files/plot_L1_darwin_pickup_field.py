
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from MITgcmutils import mds
import cmocean.cm as cm
import argparse
import sys

def read_darwin_pickup_file(input_dir,pickup_iteration):

    global_data, _, global_metadata = mds.rdmds(os.path.join(input_dir, 'pickup_darwin.'+'{:010d}'.format(pickup_iteration)), returnmeta=True)
    print(' Reading from '+os.path.join(input_dir, 'pickup_darwin.'+'{:010d}'.format(pickup_iteration)))

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)

def create_darwin_pickup_plot(config_dir, L1_model_name, pickup_iteration):

    rotate = False

    input_dir = os.path.join(config_dir,'L1', L1_model_name, 'input')

    var_names, row_bounds, var_grids, global_metadata = read_darwin_pickup_file(input_dir,pickup_iteration)


    fig = plt.figure(figsize=(8,6))
    plt.style.use("dark_background")

    counter = 1

    ff = 0

    print(' Plotting ' + var_names[ff])
    var_grid_subset = var_grids[ff][0,:,:]

    cmap = 'turbo'

    vmin = np.min(var_grid_subset[var_grid_subset != 0])
    vmax = np.max(var_grid_subset[var_grid_subset != 0])
    if vmin == vmax:
        vmin = -0.1
        vmax = 0.1

    C = plt.imshow(var_grid_subset, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    plt.colorbar(C, fraction=0.04, pad=0.04)
    plt.title(var_names[ff])

    counter += 1

    output_file = os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files',
                               L1_model_name + '_darwin_pickup_field.png')
    plt.savefig(output_file, bbox_inches='tight')
    plt.close(fig)


########################################################################################################################
# User inputs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)

