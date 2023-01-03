
import os
import argparse
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt


def read_L3_bathy_from_grid(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    depth = ds.variables['Depth'][:, :]
    ds.close()
    return(depth)

def write_exch2_files(config_dir, level_name, model_name, rows, cols,  blank_list):
    exch2_output = ' &W2_EXCH2_PARM01\n'
    exch2_output += '  W2_mapIO   = 1,\n'
    exch2_output += '  preDefTopol = 1,\n'
    # exch2_output += '  dimsFacets = ' + str(np.shape(bathy_faces[1])[1]) + ', '
    # exch2_output += str(np.shape(bathy_faces[1])[0]) + ', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n'
    exch2_output += '  dimsFacets = '+str(cols)+', '
    exch2_output += str(rows) + ',\n'
    exch2_output += '#\n'
    exch2_output += '  blankList ='
    for counter in blank_list:
        exch2_output += '  '+str(counter)+',\n'
    exch2_output += ' &'

    output_file = os.path.join(config_dir,level_name, model_name,'namelist','data.exch2')
    f = open(output_file,'w')
    f.write(exch2_output)
    f.close()

def create_exch2_file(config_dir, level_name, model_name, sNx, sNy):

    print('Creating the data.exch2 file for the '+model_name+' model')

    depth = read_L3_bathy_from_grid(config_dir, model_name)
    rows = np.shape(depth)[0]
    cols = np.shape(depth)[1]

    n_tile_rows = int(np.shape(depth)[0] / sNy)
    n_tile_cols = int(np.shape(depth)[1] / sNx)
    n_blank_cells = 0
    total_cells = 0
    blank_list = []
    for j in range(n_tile_rows):
        for i in range(n_tile_cols):
            total_cells += 1
            tile_subset = depth[j * sNy:(j + 1) * sNy, i * sNx:(i + 1) * sNx]
            if not np.any(tile_subset != 0):
                n_blank_cells += 1
                blank_list.append(total_cells)

    print('    - Domain size: '+str(np.shape(depth)))
    print('    - Tile size: '+str(sNx)+' by '+str(sNy))
    print('    - Tile rows / cols: '+str(np.shape(depth)[0] / sNy)+' / '+str(np.shape(depth)[1] / sNx))
    print('    - Blank cells: '+str(n_blank_cells))
    print('    - Non-blank cells: '+str(total_cells-n_blank_cells))
    print('    - Total cells: '+str(total_cells))
    print('    - Proc breakdown: ')
    print('        - Ivy procs: ' + str((total_cells - n_blank_cells) / 20))
    print('        - Has procs: ' + str((total_cells - n_blank_cells) / 24))
    print('        - Bro procs: ' + str((total_cells - n_blank_cells) / 28))
    print('    - Blank list: '+str(blank_list))

    write_exch2_files(config_dir, level_name, model_name, rows, cols, blank_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-l", "--level_name", action="store",
                        help="The level of the model (e.g. L1).", dest="level_name",
                        type=str, required=True)

    parser.add_argument("-m", "--model_name", action="store",
                        help="The name of the model (e.g. L1_East_Pacific).", dest="model_name",
                        type=str, required=True)

    parser.add_argument("-sNx", "--sNx", action="store",
                        help="The size of the processing tile in the x direction.", dest="sNx",
                        type=int, required=True)

    parser.add_argument("-sNy", "--sNy", action="store",
                        help="The size of the processing tile in the y direction.", dest="sNy",
                        type=int, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    model_name = args.model_name
    level_name = args.level_name
    sNx = args.sNx
    sNy = args.sNy

    create_exch2_file(config_dir, level_name, model_name, sNx, sNy)