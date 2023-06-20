import os
import simplegrid as sg
import numpy as np
from scipy.interpolate import griddata
#from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
import matplotlib.pyplot as plt
import matplotlib.path as mplP
import argparse
import ast

def map_the_bathymetry_to_subdomain(bathy_grid_faces,XC_faces,YC_faces,XC_subset,YC_subset,llc):

    bathy_45 = np.concatenate((bathy_grid_faces[4].T, bathy_grid_faces[5].T), axis=1)
    XC_45 = np.concatenate((XC_faces[4].T, XC_faces[5].T), axis=1)
    YC_45 = np.concatenate((YC_faces[4].T, YC_faces[5].T), axis=1)

    points = np.hstack([np.reshape(XC_45,(np.size(XC_45),1)),
                        np.reshape(YC_45,(np.size(YC_45),1))])
    values = np.reshape(bathy_45, (np.size(bathy_45), 1))

    grid = griddata(points,values,(XC_subset,YC_subset),method='nearest')

    return(grid)


def create_bathymetry_file(config_dir, ecco_dir):

    print('Creating the L1_1080 bathymetry file')
    llc = 1080

    # f = open(os.path.join(config_dir, 'domain_sizes.txt'))
    # dict_str = f.read()
    # f.close()
    # size_dict = ast.literal_eval(dict_str)
    # L1_size = size_dict['L1_'+str(llc)]
    # n_rows = L1_size[0]
    # n_cols = L1_size[1]
    n_rows = 360
    n_cols=480

    print('    - Reading in the LLC1080 grid tiles')
    # read in the LLC grid subset to faces
    grid_file_dir = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles')
    XC_faces = {}
    YC_faces = {}
    for i in range(1, 6):
        if i < 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc,
                                                   3 * llc)
        if i == 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc,
                                                   llc)
        if i > 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), 3 * llc,
                                                   llc)
        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face

    print('    - Reading in the LLC1080 bathymetry')
    # read in the LLC bathy subset to faces
    bathy_file = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files', 'input_init', 'bathy_llc'+str(llc))
    bathy_grid = np.fromfile(bathy_file,'>f4')
    bathy_grid_faces = {}
    points_counted = 0
    for i in range(1, 6):
        if i < 3:
            grid_subset = bathy_grid[points_counted:points_counted+3*llc*llc]
            bathy_grid_faces[i] = np.reshape(grid_subset,(3*llc,llc))
            points_counted += 3*llc*llc
        if i == 3 or i==6:
            grid_subset = bathy_grid[points_counted:points_counted + llc * llc]
            bathy_grid_faces[i] = np.reshape(grid_subset,(llc, llc))
            points_counted += llc * llc
        if i > 3 and i < 6:
            grid_subset = bathy_grid[points_counted:points_counted + 3 * llc * llc]
            bathy_grid_faces[i] = np.reshape(grid_subset,(llc, llc *3))
            points_counted += 3 * llc * llc

    print('    - Reading in the L1_1080 grid tile')
    # read in the grid subset to interpolate onto
    mitgrid_file = os.path.join(config_dir,'L1','L1_East_Pacific','input','L1_East_Pacific.mitgrid')
    grid_dict = sg.gridio.read_mitgridfile(mitgrid_file,n_cols,n_rows)
    XC_subset = grid_dict['XC'].T
    YC_subset = grid_dict['YC'].T

    # # map the LLC files into the subdomain
    print('    - Mapping the LLC'+str(llc)+' bathymetry onto the L1_1080 domain')
    bathy_grid = map_the_bathymetry_to_subdomain(bathy_grid_faces,XC_faces,YC_faces,XC_subset,YC_subset,llc)
    bathy_grid = bathy_grid[:,:,0]

    print('    - Manually filling in land cells below sea level which will cause rStar issues')
    filled_bathy_grid = np.copy(bathy_grid)
    fig = plt.figure()
    cs = plt.contour(bathy_grid, levels=[-0.5])
    plt.close(fig)
    paths = cs.collections[0].get_paths()
    path = paths[0]
    v = path.vertices
    x = v[:, 0]
    y = v[:, 1]
    x = np.concatenate([x, np.array([480]), np.array([480]), np.array([x[0]])])
    y = np.concatenate([y, np.array([0]), np.array([360]), np.array([360])])
    long_x = np.arange(480)
    long_y = np.arange(360)
    land_bound = np.column_stack([x, y])
    mplPath = mplP.Path(land_bound)
    for i in long_x:
        for j in long_y:
            if bathy_grid[j, i] > -10 and bathy_grid[j, i] < 10:
                # plt.plot(i,j,'k.')
                if mplPath.contains_point([i, j]):
                    filled_bathy_grid[j, i] = 0
    filled_bathy_grid[224:260, 326:340] = 0
    filled_bathy_grid[260:270, 273:283] = 0
    filled_bathy_grid[179:203, 335:366] = 0
    filled_bathy_grid[264, 268] = 0
    filled_bathy_grid[264, 269] = 0

    print('    - Storing as L1_1080/input/bathymetry.bin')
    output_dir = os.path.join(config_dir,'L1','L1_East_Pacific', 'input')
    output_file = os.path.join(output_dir, 'L1_East_Pacific_bathymetry.bin')
    filled_bathy_grid.ravel('C').astype('>f4').tofile(output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_path = args.ecco_path

    create_bathymetry_file(config_dir,ecco_path)
