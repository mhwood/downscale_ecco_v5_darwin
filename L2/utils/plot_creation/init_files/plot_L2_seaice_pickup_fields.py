
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from MITgcmutils import mds
import cmocean.cm as cm
import argparse


def read_seaice_pickup_file(input_dir,pickup_iteration):

    global_data, _, global_metadata = mds.rdmds(os.path.join(input_dir, 'pickup_seaice.'+'{:010d}'.format(pickup_iteration)), returnmeta=True)
    print(' Reading from '+os.path.join(input_dir, 'pickup_seaice.'+'{:010d}'.format(pickup_iteration)))

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


def create_seaice_pickup_plot(config_dir, L2_model_name, pickup_iteration):

    input_dir = os.path.join(config_dir,'L2', L2_model_name, 'input')
    output_dir = os.path.join(config_dir,'L2', L2_model_name, 'plots', 'init_files')

    var_names, row_bounds, var_grids, global_metadata = read_seaice_pickup_file(input_dir,pickup_iteration)

    fig = plt.figure(figsize=(12, 8))
    plt.style.use("dark_background")

    counter = 1

    for ff in range(len(var_names)):

        plt.subplot(2, 3, counter)

        var_grid_subset = var_grids[ff][0, :, :]

        cmap = 'viridis'

        if var_names[ff] in ['siUICE','siVICE']:
            cmap = cm.balance
            vmin = -0.15
            vmax = 0.15
        else:
            if np.any(var_grids[ff][0, :, :])>0:
                vmin = np.min(var_grid_subset[var_grids[ff][0, :, :] != 0])
                vmax = np.max(var_grid_subset[var_grids[ff][0, :, :] != 0])
            else:
                vmin=-0.1
                vmax=0.1

        if var_names[ff] in ['siAREA','siHEFF','siHSNOW']:
            cmap = cm.ice
        if var_names[ff] in ['siTICE']:
            cmap = cm.thermal

        C = plt.imshow(var_grid_subset, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        plt.colorbar(C)
        plt.title(var_names[ff])
        counter += 1

    plt.savefig(os.path.join(output_dir, L2_model_name+'_seaice_pickup_fields.png'), bbox_inches='tight')
    plt.close(fig)


########################################################################################################################
# User inputs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)

