
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from MITgcmutils import mds
import cmocean.cm as cm
import argparse
import sys

def read_L1_grid_geometry(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    AngleSN = ds.variables['AngleSN'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    ds.close()
    return(AngleCS, AngleSN)

def rotate_oriented_grids_to_natural_grids(var_names, var_grids, AngleCS, AngleSN):

    def rotate_velocity_vectors_to_natural(angle_cos, angle_sin, uvel, vvel):
        zonal_velocity = np.zeros_like(uvel)
        meridional_velocity = np.zeros_like(vvel)
        zonal_velocity[:,:] = angle_cos * uvel[:,:] - angle_sin * vvel[:,:]
        meridional_velocity[:,:] = angle_sin * uvel[:,:] + angle_cos * vvel[:,:]
        return (zonal_velocity, meridional_velocity)

    uvel_grid = var_grids[var_names.index('siUICE')]
    vvel_grid = var_grids[var_names.index('siVICE')]
    natural_uvel_grid, natural_vvel_grid = rotate_velocity_vectors_to_natural(AngleCS, AngleSN, uvel_grid, vvel_grid)
    var_grids[var_names.index('siUICE')] = natural_uvel_grid
    var_grids[var_names.index('siVICE')] = natural_vvel_grid

    return(var_grids)

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

def create_seaice_pickup_plot(config_dir, L1_model_name, pickup_iteration):

    rotate = False

    input_dir = os.path.join(config_dir,'L1', L1_model_name, 'input')

    AngleCS, AngleSN = read_L1_grid_geometry(config_dir, L1_model_name)

    var_names, row_bounds, var_grids, global_metadata = read_seaice_pickup_file(input_dir,pickup_iteration)

    if rotate:
        stitched_var_grids = rotate_oriented_grids_to_natural_grids(var_names, stitched_var_grids, AngleCS, AngleSN)

    fig = plt.figure(figsize=(12, 8))
    plt.style.use("dark_background")

    counter = 1

    for ff in range(len(var_names)):

        plt.subplot(2, 3, counter)

        print(' Plotting ' + var_names[ff])
        var_grid_subset = var_grids[ff][0,:,:]

        cmap = cm.ice

        if var_names[ff] in ['siUICE','siVICE']:
            cmap = cm.balance
            u_subset = var_grids[var_names.index('siUICE')][0,:,:]
            v_subset = var_grids[var_names.index('siVICE')][0, :, :]
            if np.any(u_subset!=0):
                vmin_u = np.min(u_subset[u_subset != 0])
                vmax_u = np.max(u_subset[u_subset != 0])
            else:
                vmax_u = 0
                vmin_u = 0
            if np.any(v_subset!=0):
                vmin_v = np.min(v_subset[v_subset != 0])
                vmax_v = np.max(v_subset[v_subset != 0])
            else:
                vmax_v = 0
                vmin_v = 0
            vmin = -np.max([np.abs(vmax_u),np.abs(vmin_u),np.abs(vmax_v),np.abs(vmin_v)])
            vmax = np.max([np.abs(vmax_u),np.abs(vmin_u),np.abs(vmax_v),np.abs(vmin_v)])
        else:
            vmin = np.min(var_grid_subset[var_grid_subset != 0])
            vmax = np.max(var_grid_subset[var_grid_subset != 0])
        if vmin == vmax:
            vmin = -0.1
            vmax = 0.1

        if var_names[ff]=='siTICE':
            cmap = cm.thermal

        C = plt.imshow(var_grid_subset, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        plt.colorbar(C, fraction=0.04, pad=0.04)
        if var_names[ff] in ['siUICE','siVICE']:
            if rotate:
                plt.title(var_names[ff]+' (Rotated)')
            else:
                plt.title(var_names[ff])
        else:
            plt.title(var_names[ff])

        counter += 1

    output_file = os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files',
                               L1_model_name + '_seaice_pickup_fields.png')
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

