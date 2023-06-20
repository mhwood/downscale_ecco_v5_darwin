
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast

def read_L3_rows_cols_from_grid(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    ds.close()
    rows = np.shape(XC)[0]
    cols = np.shape(XC)[1]
    return(rows, cols)

def create_mitgrid_plot(config_dir, L3_model_name):

    n_rows, n_cols = read_L3_rows_cols_from_grid(config_dir, L3_model_name)

    grid_keys = ['XC', 'YC', 'DXF', 'DYF', 'RAC', 'XG', 'YG', 'DXV', 'DYU',
                 'RAZ', 'DXC', 'DYC', 'RAW', 'RAS', 'DXG', 'DYG']

    file_path = os.path.join(config_dir, 'mitgrids', L3_model_name+'.mitgrid')
    entire_grid = np.fromfile(file_path, dtype='>f8')
    entire_grid = np.reshape(entire_grid, (16, n_rows + 1, n_cols + 1))

    fig = plt.figure(figsize=(12, 12))

    plt.style.use('dark_background')

    for i in range(16):
        subset_grid = entire_grid[i,:,:]
        field_name = grid_keys[i]

        if field_name in ['XC', 'YC', 'DXF', 'DYF', 'RAC']:
            subset_grid = subset_grid[:-1,:-1]
        if field_name in ['DXV']:
            subset_grid = subset_grid[:,1:-1]
        if field_name in ['DYU']:
            subset_grid = subset_grid[1:-1,:]
        if field_name in ['RAZ']:
            subset_grid = subset_grid[1:-1,1:-1]
        if field_name in ['DXC','RAW']:
            subset_grid = subset_grid[:-1,1:-1]
        if field_name in ['DYC','RAS']:
            subset_grid = subset_grid[1:-1,:-1]
        if field_name in ['DXG']:
            subset_grid = subset_grid[:,:-1]
        if field_name in ['DYG']:
            subset_grid = subset_grid[:-1,:]

        # vmin = np.min(entire_grid)
        # vmax = np.max(entire_grid)

        plt.subplot(4, 4, i+1)
        C = plt.imshow(subset_grid,origin='lower')#,vmin=vmin,vmax=vmax)
        plt.colorbar(C)
        plt.title(grid_keys[i])

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

    plt.suptitle('mitgrid file components: '+field_name+' ('+str(n_rows)+' rows by '+str(n_cols)+' columns)')

    output_file = os.path.join(config_dir, 'L3', L3_model_name, 'plots','init_files', L3_model_name+'mitgrid_components.png')
    plt.savefig(output_file,bbox_inches = 'tight')
    plt.close(fig)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)


