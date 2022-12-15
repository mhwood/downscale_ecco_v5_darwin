

import os
import numpy as np
import argparse
import ast


########################################################################################################################

def split_mitgrid_file(config_dir):
    L1_model_name = 'L1_mac_delta'

    # read in the size of the domain from the SIZE.h file
    size_file = os.path.join(config_dir, 'L1', L1_model_name, 'code_for_grid', 'SIZE.h')
    f= open(size_file)
    lines=f.read()
    f.close()
    lines=lines.split('\n')
    for line in lines:
        line = line.split()
        if len(line)>2:
            if line[0].strip()=='&' and line[1].strip() == 'sNx':
                cols = int(line[3].strip()[:-1])
            if line[0].strip()=='&' and line[1].strip() == 'sNy':
                rows = int(line[3].strip()[:-1])
            if line[0].strip()=='&' and line[1].strip() == 'nPx':
                px = int(line[3].strip()[:-1])
            if line[0].strip()=='&' and line[1].strip() == 'nPy':
                py = int(line[3].strip()[:-1])

    n_rows = rows * py
    n_cols = cols * px

    print('    - Read the following data from the SIZE.h file:')
    print('        - sNx: ' + str(cols))
    print('        - sNy: ' + str(rows))
    print('        - nPx: ' + str(px))
    print('        - nPy: ' + str(py))
    print('        - rows: ' + str(n_rows))
    print('        - cols: ' + str(n_cols))

    mitgrid_file = os.path.join(config_dir,'L1',L1_model_name,'input',L1_model_name+'.mitgrid')
    entire_grid = np.fromfile(mitgrid_file,'>f8')
    entire_grid = np.reshape(entire_grid,(16,n_rows+1,n_cols+1))

    nr = int(n_rows/rows)
    nc = int(n_cols/cols)
    counter = 1
    for ri in range(nr):
        for ci in range(nc):
            tile_subset = entire_grid[:,rows*ri:rows*(ri+1)+1,
                                        cols*ci:cols*(ci+1)+1]
            output_file = os.path.join(config_dir,'L1',L1_model_name,'input','tile'+'{:03d}'.format(counter)+'.mitgrid')
            print('     '+output_file+' shape: '+str(np.shape(tile_subset)))
            tile_subset.ravel('C').astype('>f8').tofile(output_file)
            counter+=1

    print('    - Total domain shape:' +str(np.shape(entire_grid)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1_1080 and L2_2160 directories are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    split_mitgrid_file(config_dir)