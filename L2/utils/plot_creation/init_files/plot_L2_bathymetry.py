
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import argparse
import cmocean.cm as cm
import ast

def read_L3_bathy_and_wetgrid_from_grid(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    depth = ds.variables['Depth'][:, :]
    wet_grid = ds.variables['HFacC'][:, :, :]
    ds.close()
    wet_grid = wet_grid[0, :, :]
    wet_grid[wet_grid>0] = 1
    return(depth, wet_grid)

def create_bathymetry_plot(config_dir, L3_model_name):

    depth, wet_grid = read_L3_bathy_and_wetgrid_from_grid(config_dir, L3_model_name)

    sNx = 53
    sNy = 30

    fig = plt.figure(figsize=(18, 6))
    plt.style.use('dark_background')

    plt.subplot(1, 2, 1)
    plt.contour(depth,levels=[250,500,750,1000],colors='k',linewidths=0.25)
    C = plt.imshow(depth,origin='lower',cmap = cm.deep)#,vmin=vmin,vmax=vmax)

    n_tile_rows = int(np.shape(depth)[0] / sNy)
    n_tile_cols = int(np.shape(depth)[1] / sNx)
    n_blank_cells = 0
    total_cells = 0
    for j in range(n_tile_rows):
        for i in range(n_tile_cols):
            total_cells += 1
            tile_subset = depth[j*sNy:(j+1)*sNy,i*sNx:(i+1)*sNx]
            if np.any(tile_subset!=0):
                rect = Rectangle((i*sNx,j*sNy),sNx,sNy,facecolor='none',edgecolor='k')
            else:
                rect = Rectangle((i * sNx, j * sNy), sNx, sNy, facecolor='none', edgecolor='k', hatch = '//')
                n_blank_cells+=1
            plt.text((i+0.5)*sNx,(j+0.5)*sNy,str(total_cells),ha='center',va='center',color='k', fontsize = 5)
            plt.gca().add_patch(rect)


    plt.colorbar(C, fraction=0.025, pad=0.04)
    plt.title('Bathymetry (250m Contours)')
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    plt.xlabel('With '+str(sNx)+' x '+str(sNy)+' tiles, '+str(n_blank_cells)+' blank + '+str(total_cells-n_blank_cells)+
               ' non-blank =  '+str(total_cells)+' total')

    plt.subplot(1, 2, 2)
    C = plt.imshow(wet_grid, origin='lower')  # ,vmin=vmin,vmax=vmax)
    plt.colorbar(C, fraction=0.046, pad=0.04)
    plt.title('Wet Grid')
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    output_file = os.path.join(config_dir, 'L3', L3_model_name, 'plots', 'init_files', L3_model_name+'_bathymetry.png')
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


