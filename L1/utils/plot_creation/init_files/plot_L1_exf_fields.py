
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from MITgcmutils import mds
from random import randint
import cmocean.cm as cm
import sys
import argparse


def read_exf_fields(exf_dir, var_name, n_rows, n_cols, year, month, day, hour):

    if var_name in ['UWIND','VWIND','USTRESS','VSTRESS']:
        file_path = os.path.join(exf_dir,var_name, 'L1_exf_' + var_name + '.'+str(year)+'{:02d}'.format(month)+'_rotated.bin')
    else:
        file_path = os.path.join(exf_dir, var_name,
                                 'L1_exf_' + var_name + '.' + str(year) + '{:02d}'.format(month) + '.bin')
    exf_field = np.fromfile(file_path, '>f4')

    if month in [1, 3, 5, 7, 8, 10, 12]:
        nDays = 31
    elif month in [4, 6, 9, 11]:
        nDays = 30
    else:
        if year % 4 == 0:
            nDays = 29
        else:
            nDays = 28
    n_timesteps = nDays * 4

    exf_field = np.reshape(exf_field, (n_timesteps, n_rows, n_cols))

    timestep = 4*(day-1) + int((hour-3)/6)
    exf_field = exf_field[timestep, :, :]

    return(exf_field)

def var_name_to_units(var_name):
    if var_name=='ATEMP':
        units = 'K'
    elif var_name=='AQH':
        units = 'kg/kg'
    elif var_name=='LWDOWN' or var_name=='SWDOWN':
        units = 'W/m$^2$'
    elif var_name=='UWIND' or var_name=='VWIND':
        units = 'm/s'
    elif var_name=='PRECIP' or var_name=='RUNOFF':
        units = 'm/s'
    elif var_name=='ATMOSCO2' or var_name=='IRONDUST':
        units = 'g/kg'
    else:
        units = ''
    return(units)

def plot_L1_exfs(config_dir, L1_model_name, n_rows, n_cols, year, month, day, hour):

    if 'plots' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name)):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'plots'))
    if 'init_files' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name,'plots')):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'plots','init_files'))

    var_names = ['ATEMP', 'AQH', 'LWDOWN', 'SWDOWN', 'UWIND', 'VWIND', 'PRECIP', 'RUNOFF', 'ATMOSCO2', 'IRONDUST']

    output_dir = os.path.join(config_dir,'L1', L1_model_name, 'plots', 'init_files')

    fig = plt.figure(figsize=(20, 12))
    plt.style.use("dark_background")

    counter = 1
    exf_dir = os.path.join(config_dir,'L1',L1_model_name,'input','exf')

    var_grids = []
    for ff in range(len(var_names)):
        var_grid_subset = read_exf_fields(exf_dir, var_names[ff], n_rows, n_cols, year, month, day, hour)
        if np.any(var_grid_subset!=0):
            print('   - The '+var_names[ff]+' grid has values in the range '+str(np.min(var_grid_subset[var_grid_subset!=0]))+
                  ' to '+str(np.max(var_grid_subset[var_grid_subset!=0])))
        else:
            print('   - The '+var_names[ff]+' grid is equivalently zero')
        var_grids.append(var_grid_subset)

    for ff in range(len(var_names)):

        plt.subplot(2, 5, counter)
        cmap = 'viridis'

        var_grid_subset = var_grids[ff]

        if var_names[ff] in ['UWIND','VWIND']:
            cmap = cm.balance
            val = np.max(np.abs(var_grid_subset[var_grid_subset != 0]))
            vmin = -val
            vmax = val
        else:
            if np.any(var_grids[ff][ :, :])>0:
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
        if var_names[ff] in ['IRONDUST']:
            cmap = cm.thermal

        C = plt.imshow(var_grid_subset, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        cbar = plt.colorbar(C, fraction=0.046, pad=0.04, orientation = 'horizontal')
        cbar.set_label(var_name_to_units(var_names[ff]))
        plt.title(var_names[ff])
        counter += 1

    plt.suptitle('Boundary Conditions on date ' + str(year) + '/' + str(month) + '/' + str(day))

    plt.savefig(os.path.join(output_dir, L1_model_name+'_exf_fields_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)+'_'+'{:02d}'.format(hour)+'00.png'), bbox_inches='tight')
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

    plot_L1_exfs(config_dir)

