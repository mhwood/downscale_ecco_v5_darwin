
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
    drF = ds.variables['drF'][:]
    AngleSN = ds.variables['AngleSN'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    ds.close()
    Nr = len(drF)
    return(AngleCS, AngleSN, Nr)

def rotate_oriented_grids_to_natural_grids(var_names, var_grids, AngleCS, AngleSN):

    def rotate_velocity_vectors_to_natural(angle_cos, angle_sin, uvel, vvel):
        zonal_velocity = np.zeros_like(uvel)
        meridional_velocity = np.zeros_like(vvel)
        zonal_velocity[:,:] = angle_cos * uvel[:,:] - angle_sin * vvel[:,:]
        meridional_velocity[:,:] = angle_sin * uvel[:,:] + angle_cos * vvel[:,:]
        return (zonal_velocity, meridional_velocity)

    uvel_grid = var_grids[var_names.index('Uvel')]
    vvel_grid = var_grids[var_names.index('Vvel')]
    natural_uvel_grid, natural_vvel_grid = rotate_velocity_vectors_to_natural(AngleCS, AngleSN, uvel_grid, vvel_grid)
    var_grids[var_names.index('Uvel')] = natural_uvel_grid
    var_grids[var_names.index('Vvel')] = natural_vvel_grid

    gunm1_grid = var_grids[var_names.index('GuNm1')]
    gvnm1_grid = var_grids[var_names.index('GvNm1')]
    natural_gunm1_grid, natural_gvnm1_grid = rotate_velocity_vectors_to_natural(AngleCS, AngleSN, gunm1_grid, gvnm1_grid)
    var_grids[var_names.index('GuNm1')] = natural_gunm1_grid
    var_grids[var_names.index('GvNm1')] = natural_gvnm1_grid

    gunm2_grid = var_grids[var_names.index('GuNm2')]
    gvnm2_grid = var_grids[var_names.index('GvNm2')]
    natural_gunm2_grid, natural_gvnm2_grid = rotate_velocity_vectors_to_natural(AngleCS, AngleSN, gunm2_grid, gvnm2_grid)
    var_grids[var_names.index('GuNm2')] = natural_gunm2_grid
    var_grids[var_names.index('GvNm2')] = natural_gvnm2_grid

    return(var_grids)

def read_pickup_file(input_dir,Nr,pickup_iteration):

    global_data, _, global_metadata = mds.rdmds(os.path.join(input_dir, 'pickup.'+'{:010d}'.format(pickup_iteration)), returnmeta=True)
    print(' Reading from '+os.path.join(input_dir, 'pickup.'+'{:010d}'.format(pickup_iteration)))

    has_Nr = {'uvel': True, 'vvel': True, 'theta': True,
              'salt': True, 'gunm1': True, 'gvnm1': True,
              'gunm2': True, 'gvnm2': True, 'etan': False,
              'detahdt': False, 'etah': False}

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if has_Nr[var_name.strip().lower()]:
            end_row = start_row + Nr
        else:
            end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)


def create_pickup_plot(config_dir, L1_model_name, pickup_iteration):

    if 'plots' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name)):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'plots'))
    if 'init_files' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name,'plots')):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'plots','init_files'))

    rotate = False

    input_dir = os.path.join(config_dir,'L1', L1_model_name, 'input')

    AngleCS, AngleSN, Nr = read_L1_grid_geometry(config_dir, L1_model_name)

    var_names, row_bounds, var_grids, global_metadata = read_pickup_file(input_dir,Nr,pickup_iteration)

    if rotate:
        var_grids = rotate_oriented_grids_to_natural_grids(var_names, var_grids, AngleCS, AngleSN)

    fig = plt.figure(figsize=(12, 8))
    plt.style.use("dark_background")

    counter = 1

    for ff in range(len(var_names)):

        plt.subplot(3, 4, counter)

        print(' Plotting ' + var_names[ff])
        var_grid_subset = var_grids[ff][0,:,:]

        cmap = cm.tempo_r

        if var_names[ff] in ['Uvel', 'Vvel', 'GuNm1', 'GuNm2', 'GvNm1', 'GvNm2']:
            cmap = cm.balance
            if var_names[ff] in ['Uvel', 'Vvel']:
                var_grid_subset_u = var_grids[var_names.index('Uvel')][0, :, :]
                var_grid_subset_v = var_grids[var_names.index('Vvel')][0, :, :]
                min_vel = np.min([np.min(var_grid_subset_u),np.min(var_grid_subset_v)])
                max_vel = np.max([np.max(var_grid_subset_u), np.max(var_grid_subset_v)])
                vmin = -np.max([np.abs(min_vel),np.abs(max_vel)])
                vmax = np.max([np.abs(min_vel),np.abs(max_vel)])
            if var_names[ff] in ['GuNm1', 'GvNm1']:
                var_grid_subset_u = var_grids[var_names.index('GuNm1')][0, :, :]
                var_grid_subset_v = var_grids[var_names.index('GvNm1')][0, :, :]
                min_vel = np.min([np.min(var_grid_subset_u), np.min(var_grid_subset_v)])
                max_vel = np.max([np.max(var_grid_subset_u), np.max(var_grid_subset_v)])
                vmin = -np.max([np.abs(min_vel), np.abs(max_vel)])
                vmax = np.max([np.abs(min_vel), np.abs(max_vel)])
            if var_names[ff] in ['GuNm2', 'GvNm2']:
                var_grid_subset_u = var_grids[var_names.index('GuNm2')][0, :, :]
                var_grid_subset_v = var_grids[var_names.index('GvNm2')][0, :, :]
                min_vel = np.min([np.min(var_grid_subset_u), np.min(var_grid_subset_v)])
                max_vel = np.max([np.max(var_grid_subset_u), np.max(var_grid_subset_v)])
                vmin = -np.max([np.abs(min_vel), np.abs(max_vel)])
                vmax = np.max([np.abs(min_vel), np.abs(max_vel)])
        else:
            if np.any(var_grid_subset != 0):
                vmin = np.min(var_grid_subset[var_grid_subset != 0])
                vmax = np.max(var_grid_subset[var_grid_subset != 0])
                if vmin == vmax:
                    vmin += -0.1
                    vmax += 0.1
            else:
                vmin = -0.1
                vmax = 0.1

        if var_names[ff]=='Theta':
            cmap = cm.thermal
        if var_names[ff]=='Salt':
            cmap = cm.haline

        C = plt.imshow(var_grid_subset, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        plt.colorbar(C, fraction=0.04, pad=0.04)
        if var_names[ff] in ['Uvel', 'Vvel', 'GuNm1', 'GuNm2', 'GvNm1', 'GvNm2']:
            if rotate:
                plt.title(var_names[ff]+' (Rotated)')
            else:
                plt.title(var_names[ff])
        else:
            plt.title(var_names[ff])

        counter += 1

    output_file = os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files',
                               L1_model_name + '_pickup_fields.png')
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

