
import os
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import cmocean.cm as cm


def read_boundary_condition(config_dir, L1_model_name, boundary, field_name, n_rows, n_cols, Nr, year, month, day):

    if field_name in ['UVEL', 'VVEL', 'UICE', 'VICE']:
        suffix = '_rotated.bin'
    else:
        suffix = '.bin'

    boundary_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',boundary,field_name,
                                 'L1_BC_'+boundary+'_' + field_name.upper() + '.'+str(year)+'{:02d}'.format(month)+suffix)

    boundary_grid = np.fromfile(boundary_file, '>f4')

    if month in [1, 3, 5, 7, 8, 10, 12]:
        nDays = 31
    elif month in [4, 6, 9, 11]:
        nDays = 30
    else:
        if year % 4 == 0:
            nDays = 29
        else:
            nDays = 28

    if boundary == 'north' or boundary == 'south':
        boundary_grid = boundary_grid.reshape((nDays, Nr, n_cols))
    if boundary == 'east' or boundary == 'west':
        boundary_grid = boundary_grid.reshape((nDays, Nr, n_rows))

    boundary_grid = boundary_grid[day-1,:,:]

    return(boundary_grid)

def read_boundary_masks(config_dir, L1_model_name, boundaries):
    C_masks = []
    W_masks = []
    S_masks = []

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',L1_model_name+'_grid.nc'))
    HFacC = ds.variables['HFacC'][:, :, :]
    HFacS = ds.variables['HFacS'][:, :, :]
    HFacW = ds.variables['HFacW'][:, :, :]
    ds.close()

    for boundary in boundaries:
        if boundary=='north':
            C_masks.append(HFacC[:,-1, :])
            S_masks.append(HFacS[:, -2, :])
            W_masks.append(HFacW[:, -1, :])
        if boundary=='west':
            C_masks.append(HFacC[:, :, 0])
            S_masks.append(HFacS[:, :, 0])
            W_masks.append(HFacW[:, :, 1])
        if boundary=='south':
            C_masks.append(HFacC[:, 0, :])
            S_masks.append(HFacS[:, 1, :])
            W_masks.append(HFacW[:, 0, :])
        if boundary=='east':
            C_masks.append(HFacC[:, :, -1])
            S_masks.append(HFacS[:, :, -1])
            W_masks.append(HFacW[:, :, -2])

    return(C_masks,W_masks,S_masks)

def tracer_number_to_metadata(tracer_number):
    tracer_dict = {1:['DIC','$\mu$M C',cm.amp_r],
                   2:['NO$_3$','$\mu$M N',cm.dense],
                   3:['NO$_2$','$\mu$M N',cm.dense],
                   4:['NH$_4$','$\mu$M N',cm.dense],
                   5:['PO$_4$','$\mu$M P',cm.matter],
                   6:['Fe$_T$','$\mu$M Fe',cm.turbid_r],
                   7:['SiO$_2$','$\mu$M Si',cm.gray],
                   8:['DOC','$\mu$M C',cm.amp_r],
                   9:['DON','$\mu$M N',cm.dense],
                   10:['DOP','$\mu$M P',cm.matter],
                   11:['DOFe','$\mu$M Fe',cm.turbid_r],
                   12:['POC','$\mu$M C',cm.amp_r],
                   13:['PON','$\mu$M N',cm.dense],
                   14:['POP','$\mu$M P',cm.matter],
                   15:['POFe','$\mu$M Fe',cm.turbid_r],
                   16:['POSi','$\mu$M Si',cm.gray],
                   17:['PIC','$\mu$M C',cm.amp_r],
                   18:['Alk','$\mu$M meq/m$^3$',cm.amp_r],
                   19:['O$_2$','$\mu$M O',cm.tempo_r],
                   20:['c01','$\mu$M C',cm.amp_r],
                   21:['c02','$\mu$M C',cm.amp_r],
                   22:['c03','$\mu$M C',cm.amp_r],
                   23:['c04','$\mu$M C',cm.amp_r],
                   24:['c05','$\mu$M C',cm.amp_r],
                   25:['c06','$\mu$M C',cm.amp_r],
                   26:['c07','$\mu$M C',cm.amp_r],
                   27:['Chl01','mg/m$^3$',cm.algae],
                   28:['Chl02','mg/m$^3$',cm.algae],
                   29:['Chl03','mg/m$^3$',cm.algae],
                   30:['Chl04','mg/m$^3$',cm.algae],
                   31:['Chl05','mg/m$^3$',cm.algae]}
    return(tracer_dict[tracer_number])

def create_2D_BC_plot(config_dir, L1_model_name, var_name,
                      boundaries, boundary_grids,year, month, day,
                      C_masks, W_masks, S_masks):

    min_val = 1e22
    max_val = -1e22
    for bn in range(len(boundaries)):
        boundary_grid = boundary_grids[bn]

        if np.min(boundary_grid[boundary_grid != 0]) < min_val:
            min_val = np.min(boundary_grid[boundary_grid != 0])
        if np.max(boundary_grid[boundary_grid != 0]) > max_val:
            max_val = np.max(boundary_grid[boundary_grid != 0])

    if min_val<max_val:
        val_range = max_val-min_val
    else:
        val_range = 0.1

    vmin = -0.5
    vmax = 0.5

    if var_name=='UICE':
        title_var_name = 'Zonal Seaice Ice Velocity'
        units = 'm/s'
    if var_name=='VICE':
        title_var_name = 'Meridional Seaice Ice Velocity'
        units = 'm/s'

    if var_name == 'AREA':
        vmin = -0.05
        vmax = 1.05
        title_var_name = 'Sea Ice Area'
        units = 'm$^2$/m$^2$'
    if var_name == 'HEFF':
        vmin = -0.05
        vmax = 3
        title_var_name = 'Sea Ice Effective Thickness'
        units = 'm'
    if var_name == 'HSNOW':
        vmin = -0.05
        vmax = 2
        title_var_name = 'Sea Ice Snow Effective Thickness'
        units = 'm'

    fig = plt.figure(figsize=(8, 10))
    plt.style.use('dark_background')

    for bn in range(len(boundaries)):

        plt.subplot(len(boundaries), 1, bn + 1)
        d = np.arange(np.shape(boundary_grids[bn])[1])

        if var_name == 'VICE':
            mask = S_masks[bn]
        elif var_name == 'UICE':
            mask = W_masks[bn]
        else:
            mask = C_masks[bn]
        mask = mask[:1, :np.shape(boundary_grids[bn])[1]]

        plot_line = boundary_grids[bn].ravel()
        plot_line[mask.ravel()==0] = np.nan
        plt.plot(d, plot_line)

        plt.gca().set_ylim([min_val-0.1*val_range, max_val+0.1*val_range])

        # mask_grid = np.ma.masked_where(mask != 0, mask)
        # plt.pcolormesh(X, D, mask_grid, cmap="Greys", vmin=-1, vmax=1, shading='nearest')

        plt.grid(linestyle='--',linewidth=0.5, alpha=0.5)

        plt.ylabel(boundaries[bn].capitalize() + ' Boundary\n('+units+')')
        if bn == 0:
            plt.title(title_var_name + ' Boundary Conditions on date ' + str(year) + '/' + str(month) + '/' + str(day))

        if bn == len(boundaries) - 1:
            plt.xlabel('Points Along Boundary')

    output_file = os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files', 'BCs',
                               L1_model_name + '_BC_' + var_name + '_' + str(year) + '{:02d}'.format(
                                   month) + '{:02d}'.format(day) + '.png')
    plt.savefig(output_file)
    plt.close(fig)

def create_3D_BC_plot(config_dir, L1_model_name, var_name,
                      boundaries, boundary_grids,year, month, day,
                      C_masks, W_masks, S_masks):

    min_val = 1e22
    max_val = -1e22
    for bn in range(len(boundaries)):
        boundary_grid = boundary_grids[bn]

        if np.min(boundary_grid[boundary_grid!=0])< min_val:
            min_val = np.min(boundary_grid[boundary_grid!=0])
        if np.max(boundary_grid[boundary_grid!=0])> max_val:
            max_val = np.max(boundary_grid[boundary_grid!=0])

    if var_name == 'THETA':
        cmap = cm.thermal
        vmin = min_val
        vmax = max_val
        title_var_name = 'Theta'
        units = '$^{\circ}C$'

    if var_name == 'SALT':
        cmap = cm.haline
        vmin = min_val
        vmax = max_val
        title_var_name = 'Salt'
        units = 'psu'

    if var_name == 'UVEL' or var_name=='VVEL':
        cmap = cm.balance
        vmin = -0.15
        vmax = 0.15
        if var_name=='UVEL':
            title_var_name = 'Uvel'
        if var_name=='VVEL':
            title_var_name = 'Vvel'
        units = 'm/s'

    if var_name[:6]=='PTRACE':
        tracer_number = int(var_name[-2:])
        metadata = tracer_number_to_metadata(tracer_number)
        cmap = metadata[2]
        title_var_name = metadata[0]
        vmin = min_val
        vmax = max_val
        units = metadata[1]

    fig = plt.figure(figsize=(8, 10))
    plt.style.use('dark_background')

    for bn in range(len(boundaries)):

        plt.subplot(len(boundaries), 1, bn+1)
        x = np.arange(np.shape(boundary_grids[bn])[1])
        d = np.arange(np.shape(boundary_grids[bn])[0])
        X, D = np.meshgrid(x,d)
        C = plt.pcolormesh(X, D, boundary_grids[bn], cmap=cmap, vmin=vmin, vmax=vmax, shading='nearest')

        if var_name=='VVEL':
            mask = S_masks[bn]
        elif var_name=='UVEL':
            mask = W_masks[bn]
        else:
            mask = C_masks[bn]
        mask = mask[:,:np.shape(boundary_grids[bn])[1]]
        mask_grid = np.ma.masked_where(mask != 0, mask)
        plt.pcolormesh(X, D, mask_grid, cmap="Greys", vmin=-1, vmax=1, shading='nearest')

        cbar = plt.colorbar(C)
        cbar.set_label(units)
        plt.gca().invert_yaxis()
        plt.ylabel(boundaries[bn].capitalize()+' Boundary\nDepth Levels')
        if bn==0:
            plt.title(title_var_name+' Boundary Conditions on date '+str(year)+'/'+str(month)+'/'+str(day))

        if bn==len(boundaries)-1:
            plt.xlabel('Points Along Boundary')

    output_file = os.path.join(config_dir,'L1',L1_model_name,'plots','init_files','BCs',
                               L1_model_name+'_BC_'+var_name+'_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)+'.png')
    plt.savefig(output_file)
    plt.close(fig)

def plot_L1_BCs(config_dir, L1_model_name, boundaries, Nr, n_rows, n_cols, year, month, day):

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'plot_creation','init_files'))

    if 'plots' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name)):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'plots'))
    if 'init_files' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name,'plots')):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'plots','init_files'))
    if 'BCs' not in os.listdir(os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files')):
        os.mkdir(os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files','BCs'))

    if L1_model_name in ['L1_W_Greenland', 'L1_mac_delta']:
        var_names = ['THETA', 'SALT', 'UVEL', 'VVEL', 'UICE', 'VICE', 'HSNOW', 'HEFF', 'AREA']
    else:
        var_names = ['THETA', 'SALT', 'UVEL', 'VVEL']
    for p in range(1,32):
        var_names.append('PTRACE'+'{:02d}'.format(p))

    C_masks, W_masks, S_masks = read_boundary_masks(config_dir, L1_model_name, boundaries)

    for var_name in var_names:
        print('    - Plotting the '+var_name +' boundary conditions')
        if var_name in ['AREA','HEFF','HSNOW','UICE','VICE']:
            levels = 1
        else:
            levels = Nr

        boundary_grids = []

        for boundary in boundaries:
            boundary_grid = read_boundary_condition(config_dir, L1_model_name, boundary, var_name,
                                                    n_rows, n_cols, levels, year, month, day)
            boundary_grids.append(boundary_grid)

        if var_name in ['AREA','HEFF','HSNOW','UICE','VICE']:
            create_2D_BC_plot(config_dir, L1_model_name, var_name,
                              boundaries, boundary_grids, year, month, day,
                              C_masks, W_masks, S_masks)
        else:
            create_3D_BC_plot(config_dir, L1_model_name, var_name,
                              boundaries, boundary_grids,year, month, day,
                              C_masks, W_masks, S_masks)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_L1_BCs(config_dir)
   

