
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from MITgcmutils import mds
from random import randint
import cmocean.cm as cm
import argparse

def read_L3_rows_cols_from_grid(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    ds.close()
    rows = np.shape(XC)[0]
    cols = np.shape(XC)[1]
    return(rows, cols)


def read_exf_fields(exf_dir, var_name, n_rows, n_cols, timestep):

    file_path = os.path.join(exf_dir,'L3_exf_'+var_name+'_1992')
    exf_field = np.fromfile(file_path,'>f4')
    total_timesteps = int(np.size(exf_field)/(n_rows*n_cols))
    print('    - total timsteps: '+str(total_timesteps))
    exf_field = np.reshape(exf_field,(total_timesteps,n_rows, n_cols))
    exf_field = exf_field[timestep,:,:]

    return(exf_field)


def create_exf_plot(config_dir, L3_model_name):

    output_dir = os.path.join(config_dir,'L3', L3_model_name, 'plots', 'init_files')

    n_rows, n_cols = read_L3_rows_cols_from_grid(config_dir, L3_model_name)
    print(n_rows, n_cols)

    var_names = ['UWIND','VWIND','ATEMP','AQH','PRECIP','SWDOWN','LWDOWN','RUNOFF']
    # var_names = ['ATEMP']

    fig = plt.figure(figsize=(20, 12))
    plt.style.use("dark_background")

    timestep = 90*4+17

    counter = 1
    exf_dir = os.path.join(config_dir,'L3',L3_model_name,'input','exf')

    var_grids = []
    for ff in range(len(var_names)):
        var_grid_subset = read_exf_fields(exf_dir, var_names[ff], n_rows, n_cols, timestep)
        var_grids.append(var_grid_subset)

    for ff in range(len(var_names)):

        var_grid_subset = var_grids[ff]

        plt.subplot(2, 4, counter)
        cmap = 'viridis'

        if var_names[ff] in ['UWIND','VWIND']:
            cmap = cm.balance
            vmin = -10
            vmax = 10
        else:
            if np.any(var_grid_subset)>0:
                vmin = np.min(var_grid_subset[var_grid_subset != 0])
                vmax = np.max(var_grid_subset[var_grid_subset != 0])
            else:
                vmin=-0.1
                vmax=0.1

        if var_names[ff] in ['AQH']:
            cmap = cm.haline
        if var_names[ff] in ['RUNOFF']:
            cmap = cm.rain
        if var_names[ff] in ['PRECIP']:
            cmap = cm.rain
        if var_names[ff] in ['ATEMP']:
            cmap = cm.thermal
        if var_names[ff] in ['SWDOWN','LWDOWN']:
            cmap = cm.thermal

        C = plt.imshow(var_grid_subset, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        plt.colorbar(C)
        plt.title(var_names[ff])
        counter += 1

    plt.savefig(os.path.join(output_dir, L3_model_name+'_exf_fields_timestep_'+str(timestep)+'.png'), bbox_inches='tight')
    plt.close(fig)


########################################################################################################################
# User inputs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)

