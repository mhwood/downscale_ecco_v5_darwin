



import os
import simplegrid as sg
import numpy as np
from scipy.interpolate import griddata
import argparse
import ast


########################################################################################################################

def split_mitgrid_file(config_dir,rows,cols,for_grid):

    print('Splitting the mitgrid file')

    n_rows = 510
    n_cols = 840
    model_name = 'L2_Santa_Barbara'

    if for_grid:
        for_grid_suffix = '_for_grid'
    else:
        for_grid_suffix = ''

    # read in the grid subset to interpolate onto
    mitgrid_file = os.path.join(config_dir,'L2',model_name,'input',model_name+'.mitgrid')
    entire_grid = np.fromfile(mitgrid_file,'>f8')
    entire_grid = np.reshape(entire_grid,(16,n_rows+1,n_cols+1))

    nr = int(n_rows/rows)
    nc = int(n_cols/cols)
    counter = 1
    for ri in range(nr):
        for ci in range(nc):
            tile_subset = entire_grid[:,rows*ri:rows*(ri+1)+1,
                                        cols*ci:cols*(ci+1)+1]
            output_file = os.path.join(config_dir,'L2',model_name,'input'+for_grid_suffix,'tile'+'{:03d}'.format(counter)+'.mitgrid')
            print('     '+output_file+' shape: '+str(np.shape(tile_subset)))
            tile_subset.ravel('C').astype('>f8').tofile(output_file)
            counter+=1

    print(np.shape(entire_grid))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-f", "--for_grid", action="store",
                        help="Whether the split grid is for the for_grid dir (1 for True, 0 for False)", dest="for_grid",
                        type=int, required=False, default=0)

    parser.add_argument("-r", "--rows", action="store",
                        help="The number of rows in the tile.", dest="rows",
                        type=int, required=True)

    parser.add_argument("-c", "--cols", action="store",
                        help="The number of rows in the tile.", dest="cols",
                        type=int, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    rows = args.rows
    cols = args.cols
    for_grid = args.for_grid

    if for_grid==1:
        for_grid = True
    else:
        for_grid = False

    split_mitgrid_file(config_dir,rows,cols,for_grid)