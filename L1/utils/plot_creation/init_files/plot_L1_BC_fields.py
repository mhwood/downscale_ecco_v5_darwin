
import os
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
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
                      boundaries, boundary_grids,year, month, day):

    # min_val = np.min([np.min(north_boundary[north_boundary!=0]),
    #                   np.min(south_boundary[south_boundary!=0]),
    #                   np.min(east_boundary[east_boundary!=0]),
    #                   np.min(west_boundary[west_boundary!=0])])
    #
    # max_val = np.max([np.max(north_boundary[north_boundary != 0]),
    #                   np.max(south_boundary[south_boundary != 0]),
    #                   np.max(east_boundary[east_boundary != 0]),
    #                   np.max(west_boundary[west_boundary != 0])])

    vmin = -0.5
    vmax = 0.5

    if var_name=='AREA':
        vmin=-0.05
        vmax=1.05
    if var_name=='HEFF':
        vmin = -0.05
        vmax = 3
    if var_name=='HSNOW':
        vmin = -0.05
        vmax = 2

    fig = plt.figure(figsize=(8, 10))
    plt.style.use('dark_background')

    plt.subplot(4, 1, 1)
    C = plt.plot(north_boundary.ravel())
    plt.grid(linestyle='--',alpha=0.5)
    plt.ylabel('North')
    plt.title(var_name+' Boundary Conditions at timestep = '+str(timestep))
    plt.gca().set_ylim([vmin,vmax])

    plt.subplot(4, 1, 2)
    C = plt.plot(south_boundary.ravel())
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('South')
    plt.gca().set_ylim([vmin, vmax])

    plt.subplot(4, 1, 3)
    C = plt.plot(west_boundary.ravel())
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('West')
    plt.gca().set_ylim([vmin, vmax])

    plt.subplot(4, 1, 4)
    C = plt.plot(east_boundary.ravel())
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('East')
    plt.xlabel('Points Along Boundary')
    plt.gca().set_ylim([vmin, vmax])

    output_file = os.path.join(config_dir,'L1_grid',L1_model_name,'plots','init_files',L1_model_name+'_BC_'+var_name+'_'+str(timestep)+'.png')
    plt.savefig(output_file)
    plt.close(fig)


def create_3D_BC_plot(config_dir, L1_model_name, var_name,
                      boundaries, boundary_grids,year, month, day):

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
        plot_grid = np.ma.masked_where(boundary_grids[bn]==0,boundary_grids[bn])
        C = plt.pcolormesh(X, D, plot_grid, cmap=cmap, vmin=vmin, vmax=vmax, shading='nearest')
        cbar = plt.colorbar(C)
        cbar.set_label(units)
        plt.gca().invert_yaxis()
        plt.ylabel(boundaries[bn].capitalize()+' Boundary\nDepth Levels')
        if bn==0:
            plt.title(title_var_name+' Boundary Conditions on date '+str(year)+'/'+str(month)+'/'+str(day))

        if bn==len(boundaries)-1:
            plt.xlabel('Points Along Boundary')

    output_file = os.path.join(config_dir,'L1',L1_model_name,'plots','init_files',
                               L1_model_name+'_BC_'+var_name+'_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)+'.png')
    plt.savefig(output_file)
    plt.close(fig)

def plot_L1_BCs(config_dir, L1_model_name, boundaries, Nr, n_rows, n_cols, year, month, day):

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'plot_creation','init_files'))

    if 'plots' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name)):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'plots'))
    if 'init_files' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name,'plots')):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'plots','init_files'))

    if L1_model_name in ['L1_W_Greenland', 'L1_mac_delta']:
        var_names = ['THETA', 'SALT', 'UVEL', 'VVEL', 'UICE', 'VICE', 'HSNOW', 'HEFF', 'AREA']
    else:
        var_names = ['THETA', 'SALT', 'UVEL', 'VVEL']
    for p in range(1,32):
        var_names.append('PTRACE'+'{:02d}'.format(p))

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
                              boundaries, boundary_grids,year, month, day)
        else:
            create_3D_BC_plot(config_dir, L1_model_name, var_name,
                              boundaries, boundary_grids,year, month, day)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_L1_BCs(config_dir)
   

