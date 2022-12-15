
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast
import sys

def create_mitgrid_plot(config_dir, L1_model_name, n_rows, n_cols):


    grid_keys = ['XC', 'YC', 'DXF', 'DYF', 'RAC', 'XG', 'YG', 'DXV', 'DYU',
                 'RAZ', 'DXC', 'DYC', 'RAW', 'RAS', 'DXG', 'DYG']
    file_path = os.path.join(config_dir,'L1', L1_model_name, 'input', L1_model_name+'.mitgrid')
    stitched_entire_grid = np.fromfile(file_path, dtype='>f8')
    stitched_entire_grid = np.reshape(stitched_entire_grid, (16, n_rows+1, n_cols+1))
    stitched_entire_grid = stitched_entire_grid[:,:-1,:-1]

    fig = plt.figure(figsize=(12, 12))

    plt.style.use('dark_background')

    for i in range(16):
        subset_grid = stitched_entire_grid[i,:,:]
        field_name = grid_keys[i]

        if field_name in ['XC', 'YC', 'DXF', 'DYF', 'RAC']:
            subset_grid = subset_grid#[:-1,:-1]
        if field_name in ['DXV']:
            subset_grid = subset_grid#[:,1:-1]
        if field_name in ['DYU']:
            subset_grid = subset_grid#[1:-1,:]
        if field_name in ['RAZ']:
            subset_grid = subset_grid#[1:-1,1:-1]
        if field_name in ['DXC','RAW']:
            subset_grid = subset_grid#[:-1,1:-1]
        if field_name in ['DYC','RAS']:
            subset_grid = subset_grid#[1:-1,:-1]
        if field_name in ['DXG']:
            subset_grid = subset_grid#[:,:-1]
        if field_name in ['DYG']:
            subset_grid = subset_grid#[:-1,:]

        # vmin = np.min(entire_grid)
        # vmax = np.max(entire_grid)

        plt.subplot(4, 4, i+1)
        C = plt.imshow(subset_grid,origin='lower')#,vmin=vmin,vmax=vmax)
        plt.colorbar(C)
        plt.title(grid_keys[i])

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

    plt.suptitle('mitgrid file components: '+field_name+' ('+str(np.shape(stitched_entire_grid)[1])+' rows by '+str(np.shape(stitched_entire_grid)[2])+' columns)')

    output_file = os.path.join(config_dir, 'L1', L1_model_name, 'plots','init_files', L1_model_name+'_mitgrid_components.png')
    plt.savefig(output_file,bbox_inches = 'tight')
    plt.close(fig)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)


