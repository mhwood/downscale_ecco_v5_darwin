
import os
import numpy as np
from scipy.interpolate import griddata
import argparse
import ast
import matplotlib.pyplot as plt

def read_domain_face_geometry_from_mitgrid(config_dir, L1_model_name, n_rows, n_cols):

    file_path = os.path.join(config_dir, 'L1', L1_model_name, 'input', L1_model_name + '.mitgrid')
    stitched_entire_grid = np.fromfile(file_path, dtype='>f8')
    stitched_entire_grid = np.reshape(stitched_entire_grid, (16, n_rows + 1, n_cols + 1))
    stitched_entire_grid = stitched_entire_grid[:, :-1, :-1]

    XC = stitched_entire_grid[0,:,:]
    YC = stitched_entire_grid[1,:,:]

    return(XC, YC)

def read_llc_face_geometry(ecco_dir, face, llc=1080):

    grid_file = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles', 'tile' + '{:03d}'.format(face) + '.mitgrid')
    grid = np.fromfile(grid_file, '>f8')

    if face==1 or face==2:
        grid = grid.reshape(16,3*llc+1,llc+1)
    if face==3:
        grid = grid.reshape(16,llc+1,llc+1)
    if face==4 or face==5:
        grid = grid.reshape(16,llc+1,3*llc+1)

    XC = grid[0,:-1,:-1]
    YC = grid[1, :-1, :-1]

    return(XC, YC)

def read_llc_bathymetry_file_to_faces(grid_file,llc=1080):
    grid = np.fromfile(grid_file, '>f4')

    grid_faces = {}

    # face 1
    face_1_grid = grid[:3 * llc * llc]
    face_1_grid = np.reshape(face_1_grid, (3 * llc, llc))
    grid_faces[1] = face_1_grid

    # face 3
    face_3_grid = grid[(6) * llc * llc:(6+1) * llc * llc]
    face_3_grid = np.reshape(face_3_grid, (llc, llc))
    grid_faces[3] = face_3_grid

    # face 4
    face_4_grid = grid[(6+1) * llc * llc:(6+1+3) * llc * llc]
    face_4_grid = np.reshape(face_4_grid, (llc, 3 * llc))
    grid_faces[4] = face_4_grid

    # face 3
    face_5_grid = grid[(6+1+3) * llc * llc:(6+1+6) * llc * llc]
    face_5_grid = np.reshape(face_5_grid, (llc, 3 * llc))
    grid_faces[5] = face_5_grid

    return(grid_faces)

def generate_connected_mask(start_row, start_col, wet_grid):

    if wet_grid[start_row,start_col]==0:
        raise ValueError(' The start row/col location is dry')

    rows = np.arange(np.shape(wet_grid)[0])
    cols = np.arange(np.shape(wet_grid)[1])
    Cols,Rows = np.meshgrid(cols,rows)

    mask_grid = 1-np.copy(wet_grid)
    mask_grid[start_row,start_col] = 2
    # in the mask, 0 means unverified
    # 1 is verified dry
    # 2 is verified wet

    # plt.imshow(mask_grid)
    # plt.show()

    is_remaining = np.logical_and(mask_grid==0,wet_grid==1)
    n_remaining = np.sum(is_remaining)
    # print(n_remaining)
    continue_iter = True
    for i in range(n_remaining):
        if continue_iter:
            # get the wet rows, cols, and their current mask values
            Wet_Rows = Rows[wet_grid == 1]
            Wet_Cols = Cols[wet_grid == 1]
            Mask_Vals = mask_grid[wet_grid == 1]

            # reduce these to the ones that havent been verified yet
            Wet_Rows = Wet_Rows[Mask_Vals == 0]
            Wet_Cols = Wet_Cols[Mask_Vals == 0]
            Mask_Vals = Mask_Vals[Mask_Vals == 0]

            if len(Mask_Vals)>0:

                # for each row/col, see if its connected to one we've verified is connected
                rows_remaining,cols_remaining = np.where(is_remaining)
                for ri in range(n_remaining):
                    row = rows_remaining[ri]
                    col = cols_remaining[ri]

                    # # this bit allows for diagonal spreading
                    # row_col_dist = ((Wet_Rows.astype(float)-row)**2 + (Wet_Cols.astype(float)-col)**2)**0.5
                    # closest_index = np.argmin(row_col_dist)
                    # if row_col_dist[closest_index]<np.sqrt(2):
                    #     var_grid[row,col] = Wet_Vals[closest_index]

                    # this bit allows for only up/dow/left/right spreading
                    if row<np.shape(wet_grid)[0]-1:
                        if mask_grid[row+1,col] == 2:
                            mask_grid[row,col] = 2
                    if row > 0:
                        if mask_grid[row - 1, col] == 2:
                            mask_grid[row,col] = 2
                    if col<np.shape(wet_grid)[1]-1:
                        if mask_grid[row,col+1] == 2:
                            mask_grid[row,col] = 2
                    if col > 0:
                        if mask_grid[row, col-1] == 2:
                            mask_grid[row,col] = 2


                is_remaining = np.logical_and(mask_grid == 0, wet_grid == 1)
                n_remaining_now = np.sum(is_remaining)

                # plt.subplot(1,2,1)
                # plt.imshow(wet_grid,cmap='Greys_r')
                # plt.subplot(1, 2, 2)
                # plt.imshow(mask_grid)
                # plt.show()

                if n_remaining_now<n_remaining:
                    n_remaining = n_remaining_now
                else:
                    n_remaining = n_remaining_now
                    continue_iter=False
            else:
                continue_iter = False

    return(mask_grid)

def fill_unconnected_areas(bathy_grid):

    # C = plt.imshow(stitched_bathy,origin='lower')
    # plt.colorbar(C)
    # plt.show()

    start_row = 300
    start_col = 300

    wet_grid = (bathy_grid<0).astype(int)
    mask_grid = generate_connected_mask(start_row, start_col, wet_grid)

    # plt.imshow(mask_grid, origin='lower')
    # plt.show()

    bathy_grid[mask_grid == 0] = 0

    return(bathy_grid)

def create_L1_W_Greenland_bathymetry(config_dir, ecco_dir, print_status = 4):

    model_name = 'L1_W_Greenland'
    n_rows = 720
    n_cols = 540
    llc = 1080

    if print_status>=1:
        print('    - Reading bathymetry from the ECCO LLC'+str(llc)+' model')

    XC, YC = read_domain_face_geometry_from_mitgrid(config_dir, model_name, n_rows, n_cols)

    bathy_file = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init', 'bathy_llc'+str(llc))
    bathy_faces = read_llc_bathymetry_file_to_faces(bathy_file,llc=llc)

    all_points_started = False

    for face in [1,3,4,5]:

        llc_face_XC, llc_face_YC = read_llc_face_geometry(ecco_dir, face)

        points = np.column_stack([llc_face_XC.ravel(),llc_face_YC.ravel()])
        values = bathy_faces[face].reshape((np.size(bathy_faces[face]),1))

        if not all_points_started:
            all_points = points
            all_values = values
            all_points_started = True
        else:
            all_points = np.vstack([all_points,points])
            all_values = np.vstack([all_values,values])

    bathy_grid = griddata(all_points,all_values,(XC, YC), method='nearest')
    bathy_grid = bathy_grid[:,:,0]

    bathy_grid_unfilled = np.copy(bathy_grid)

    if print_status>=1:
        print('    - Filling in unconnected areas')

    bathy_grid_filled = np.copy(bathy_grid_unfilled)

    # # close off gulf of st lawrence
    # bathy_grid_filled[300:310, 480:510] = 0
    # # close off hudson bay
    # bathy_grid_filled[150:160, 225:280] = 0

    # close off gulf of st lawrence
    bathy_grid_filled[30:60, 300:310] = 0
    # close off hudson bay
    bathy_grid_filled[260:315, 150:160] = 0

    for i in range(434,445):
        bathy_grid_filled[i, -2] = 0
        bathy_grid_filled[i, -1] = 0

    bathy_grid_filled = fill_unconnected_areas(bathy_grid_filled)

    # # fill in some annoying edge issues
    # bathy_grid_filled[440:444, 538:] = 0

    plt.imshow(bathy_grid_filled[430:450, 500:], origin='lower', vmin=-20, vmax=0)
    plt.show()

    # plt.subplot(2,2,1)
    # C = plt.imshow(XC,origin='lower')
    # plt.colorbar(C)
    # plt.title('XC')
    #
    # plt.subplot(2, 2, 2)
    # C = plt.imshow(YC, origin='lower')
    # plt.colorbar(C)
    # plt.title('YC')
    #
    # plt.subplot(2, 2, 3)
    # C = plt.imshow(bathy_grid_unfilled, origin='lower')
    # plt.colorbar(C)
    # plt.title('Bathy (Unfilled)')
    #
    # plt.subplot(2, 2, 4)
    # C = plt.imshow(bathy_grid_filled-bathy_grid_unfilled, origin='lower')
    # plt.colorbar(C)
    # plt.title('Bathy (Filled)')
    # plt.show()

    if print_status>=1:
        print('    - Outputting bathymetry to '+model_name+'_bathymetry.bin')
    output_file = os.path.join(config_dir, 'L1', model_name, 'input', model_name+'_bathymetry.bin')
    bathy_grid_filled.ravel(order='C').astype('>f4').tofile(output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_dir", action="store",
                        help="The directory where the ECCO files are stored.", dest="ecco_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_dir = args.ecco_dir

    create_L1_W_Greenland_bathymetry(config_dir, ecco_dir)
   

