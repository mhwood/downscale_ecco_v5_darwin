
import os
import numpy as np
import simplegrid as sg
import matplotlib.pyplot as plt
import netCDF4 as nc4
from MITgcmutils import mds

def get_ecco_metadata(metadata_name):
    if metadata_name=='Nr':
        metadata = 50
    elif metadata_name == 'sNx':
        metadata = 90
    elif metadata_name == 'sNy':
        metadata = 90
    elif metadata_name == 'delR':
        metadata = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                              10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                              31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                              93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                              139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                              341.50, 364.50, 387.50, 410.50, 433.50, 456.50])
    else:
        raise ValueError('metadata_name not recognized')
    return(metadata)


#################################################################################################
# grid functions

def read_ecco_field_to_faces(file_path, llc, dim):

    Nr = 50

    grid_array = np.fromfile(file_path,'>f4')
    N = 13*llc*llc

    field_faces = {}

    if dim == 2:
        points_counted = 0
        for i in range(1, 6):
            if i < 3:
                n_points = 3 * llc * llc
                grid = grid_array[points_counted:points_counted + n_points]
                grid = np.reshape(grid, (3 * llc, llc))
            if i == 3:
                n_points = llc * llc
                grid = grid_array[points_counted:points_counted + n_points]
                grid = np.reshape(grid, (llc, llc))
            if i > 3:
                n_points = 3 * llc * llc
                grid = grid_array[points_counted:points_counted + n_points]
                grid = np.reshape(grid, (llc, 3 * llc))
            field_faces[i] = grid
            points_counted += n_points

    if dim==3:

        for i in range(1, 6):
            if i < 3:
                face_grid = np.zeros((Nr, 3 * llc, llc))
            elif i == 3:
                face_grid = np.zeros((Nr, llc, llc))
            if i > 3:
                face_grid = np.zeros((Nr, llc, 3 * llc))
            field_faces[i]=face_grid

        for nr in range(Nr):
            points_counted = 0
            level_grid = grid_array[nr * N:(nr + 1) * N]
            for i in range(1,6):
                if i < 3:
                    n_points = 3*llc*llc
                    grid = level_grid[points_counted:points_counted+n_points]
                    grid = np.reshape(grid,(3*llc,llc))
                if i == 3:
                    n_points = llc * llc
                    grid = level_grid[points_counted:points_counted + n_points]
                    grid = np.reshape(grid, (llc, llc))
                if i > 3:
                    n_points = 3 * llc * llc
                    grid = level_grid[points_counted:points_counted + n_points]
                    grid = np.reshape(grid, (llc, 3*llc))
                field_faces[i][nr,:,:] = grid

                points_counted += n_points

        # plt.imshow(grid,origin='lower')
        # plt.show()

    return(field_faces)

def read_ecco_faces_to_tiles(ecco_faces, llc, dim):
    ecco_tiles = {}
    if dim==2:
        ecco_tiles[1] = ecco_faces[1][:llc,:]
        ecco_tiles[2] = ecco_faces[1][llc:2*llc, :]
        ecco_tiles[3] = ecco_faces[1][2*llc:, :]
        ecco_tiles[4] = ecco_faces[2][:llc, :]
        ecco_tiles[5] = ecco_faces[2][llc:2*llc, :]
        ecco_tiles[6] = ecco_faces[2][2*llc:, :]
        ecco_tiles[7] = ecco_faces[3][:,:]
        ecco_tiles[8] = ecco_faces[4][:, :llc]
        ecco_tiles[9] = ecco_faces[4][:, llc:2*llc]
        ecco_tiles[10] = ecco_faces[4][:, 2*llc:]
        ecco_tiles[11] = ecco_faces[5][:, :llc]
        ecco_tiles[12] = ecco_faces[5][:, llc:2*llc]
        ecco_tiles[13] = ecco_faces[5][:, 2*llc:]
    if dim==3:
        ecco_tiles[1] = ecco_faces[1][:, :llc,:]
        ecco_tiles[2] = ecco_faces[1][:, llc:2*llc, :]
        ecco_tiles[3] = ecco_faces[1][:, 2*llc:, :]
        ecco_tiles[4] = ecco_faces[2][:, :llc, :]
        ecco_tiles[5] = ecco_faces[2][:, llc:2*llc, :]
        ecco_tiles[6] = ecco_faces[2][:, 2*llc:, :]
        ecco_tiles[7] = ecco_faces[3]
        ecco_tiles[8] = ecco_faces[4][:, :, :llc]
        ecco_tiles[9] = ecco_faces[4][:, :, llc:2*llc]
        ecco_tiles[10] = ecco_faces[4][:, :, 2*llc:]
        ecco_tiles[11] = ecco_faces[5][:, :, :llc]
        ecco_tiles[12] = ecco_faces[5][:, :, llc:2*llc]
        ecco_tiles[13] = ecco_faces[5][:, :, 2*llc:]
    return(ecco_tiles)

def read_ecco_geometry_to_faces(ecco_dir,llc):
    grid_file_dir = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles')
    XC_faces = {}
    YC_faces = {}
    for i in [1, 2, 3, 4, 5]:
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

    angleCS_path = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'input_init','AngleCS.data')
    AngleCS_faces = read_ecco_field_to_faces(angleCS_path, llc, dim=2)

    angleSN_path = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'input_init', 'AngleSN.data')
    AngleSN_faces = read_ecco_field_to_faces(angleSN_path, llc, dim=2)

    hFacC_path = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'input_init', 'hFacC.data')
    hFacC_faces = read_ecco_field_to_faces(hFacC_path, llc, dim=3)

    # plt.imshow(hFacC_faces[5][5,:,:])
    # plt.show()

    return(XC_faces, YC_faces, AngleCS_faces, AngleSN_faces, hFacC_faces)

def ecco_tile_face_row_col_bounds(tile_number, llc, sNx, sNy):

    n_faces_tiles = 27
    if tile_number <= 27:
        face = 1
    elif tile_number > 27 and tile_number <= 54:
        face = 2
    elif tile_number > 54 and tile_number <= 63:
        face = 3
        n_faces_tiles = 9
    elif tile_number > 63 and tile_number <= 90:
        face = 4
    else:
        face = 5

    face_to_first_tile = {1:1, 2:28, 3:55, 4:64, 5:91}
    face_to_dims = {1:(3*llc,llc), 2:(3*llc,llc), 3:(llc,llc), 4:(llc,3*llc), 5:(llc,3*llc)}

    face_shape = face_to_dims[face]
    first_face_tile_number = face_to_first_tile[face]

    row_counter = 0
    col_counter = 0

    for counter in range(first_face_tile_number,first_face_tile_number+n_faces_tiles):

        # print(counter, row_counter, col_counter)

        if counter == tile_number:
            tile_row = row_counter
            tile_col = col_counter
        else:
            col_counter += 1

        if col_counter==face_shape[1]/sNx:
            col_counter = 0
            row_counter += 1

    min_row = tile_row * sNy
    min_col = tile_col * sNx

    return(face, min_row, min_col)

def read_ecco_grid_geometry(ecco_dir,llc,ordered_ecco_tiles,ordered_ecco_tile_rotations):

    ecco_Nr = 50
    ecco_sNx = 90
    ecco_sNy = 90

    XC_faces, YC_faces, AngleCS_faces, AngleSN_faces, hFacC_faces = read_ecco_geometry_to_faces(ecco_dir, llc)

    ecco_XC = np.zeros((ecco_sNy*len(ordered_ecco_tiles),ecco_sNx*len(ordered_ecco_tiles[0])))
    ecco_YC = np.zeros((ecco_sNy * len(ordered_ecco_tiles), ecco_sNx * len(ordered_ecco_tiles[0])))
    ecco_AngleCS = np.zeros((ecco_sNy * len(ordered_ecco_tiles), ecco_sNx * len(ordered_ecco_tiles[0])))
    ecco_AngleSN = np.zeros((ecco_sNy * len(ordered_ecco_tiles), ecco_sNx * len(ordered_ecco_tiles[0])))
    ecco_hfacC = np.zeros((ecco_Nr, ecco_sNy * len(ordered_ecco_tiles), ecco_sNx * len(ordered_ecco_tiles[0])))

    for r in range(len(ordered_ecco_tiles)):
        for c in range(len(ordered_ecco_tiles[r])):
            ecco_tile_number = ordered_ecco_tiles[r][c]
            face, min_row, min_col = ecco_tile_face_row_col_bounds(ecco_tile_number, llc, ecco_sNx,
                                                                                     ecco_sNy)

            # print('  - Reading tile '+str(ecco_tile_number)+' from face '+str(face)+' (rows '+str(min_row)+' to '+str(min_row+ecco_sNy) + \
            #       ', cols '+str(min_col)+', '+str(min_col+ecco_sNx)+')')

            # print(ecco_tile_number,face,min_row,min_col)
            XC = XC_faces[face][min_row:min_row + ecco_sNy, min_col:min_col + ecco_sNx]
            YC = YC_faces[face][min_row:min_row + ecco_sNy, min_col:min_col + ecco_sNx]
            AngleCS = AngleCS_faces[face][min_row:min_row + ecco_sNy, min_col:min_col + ecco_sNx]
            AngleSN = AngleSN_faces[face][min_row:min_row + ecco_sNy, min_col:min_col + ecco_sNx]
            hFacC = hFacC_faces[face][:, min_row:min_row + ecco_sNy, min_col:min_col + ecco_sNx]

            # rotate things as necessary
            for n in range(ordered_ecco_tile_rotations[r][c]):
                XC = np.rot90(XC)
                YC = np.rot90(YC)
                hFacC = np.rot90(hFacC, axes=(1, 2))
                AngleCS = np.rot90(AngleCS)
                AngleSN = np.rot90(AngleSN)

           # put it into the big grid
            ecco_hfacC[:, r * ecco_sNy:(r + 1) * ecco_sNy, c * ecco_sNx:(c + 1) * ecco_sNx] = hFacC
            ecco_XC[r * ecco_sNy:(r + 1) * ecco_sNy, c * ecco_sNx:(c + 1) * ecco_sNx] = XC
            ecco_YC[r * ecco_sNy:(r + 1) * ecco_sNy, c * ecco_sNx:(c + 1) * ecco_sNx] = YC
            ecco_AngleCS[r * ecco_sNy:(r + 1) * ecco_sNy, c * ecco_sNx:(c + 1) * ecco_sNx] = AngleCS
            ecco_AngleSN[r * ecco_sNy:(r + 1) * ecco_sNy, c * ecco_sNx:(c + 1) * ecco_sNx] = AngleSN

    ecco_delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                        10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                        31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                        93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                        139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                        341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    return(ecco_XC,ecco_YC,ecco_AngleCS,ecco_AngleSN,ecco_hfacC,ecco_delR)


#################################################################################################
# exf functions

def read_ecco_runoff_file(ecco_dir, file_name, ordered_ecco_tiles, ordered_ecco_tile_rotations, llc):

    print('      - Reading from '+file_name)

    timesteps = 12
    file_path = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'era_xx', file_name)
    grid_array = np.fromfile(file_path, '>f4')
    N = 13 * llc * llc

    runoff_faces = {}

    for i in range(1, 6):
        if i < 3:
            face_grid = np.zeros((timesteps, 3 * llc, llc))
        elif i == 3:
            face_grid = np.zeros((timesteps, llc, llc))
        if i > 3:
            face_grid = np.zeros((timesteps, llc, 3 * llc))
        runoff_faces[i] = face_grid

    for nr in range(timesteps):
        points_counted = 0
        level_grid = grid_array[nr * N:(nr + 1) * N]
        for i in range(1, 6):
            if i < 3:
                n_points = 3 * llc * llc
                grid = level_grid[points_counted:points_counted + n_points]
                grid = np.reshape(grid, (3 * llc, llc))
            if i == 3:
                n_points = llc * llc
                grid = level_grid[points_counted:points_counted + n_points]
                grid = np.reshape(grid, (llc, llc))
            if i > 3:
                n_points = 3 * llc * llc
                grid = level_grid[points_counted:points_counted + n_points]
                grid = np.reshape(grid, (llc, 3 * llc))
            runoff_faces[i][nr, :, :] = grid

            points_counted += n_points

    # for face in range(1,6):
    #     plt.imshow(field_faces[face][0,:,:],origin='lower')
    #     plt.show()

    ecco_sNx = 90
    ecco_sNy = 90

    ecco_Runoff = np.zeros((timesteps, ecco_sNy * len(ordered_ecco_tiles), ecco_sNx * len(ordered_ecco_tiles[0])))

    for r in range(len(ordered_ecco_tiles)):
        for c in range(len(ordered_ecco_tiles[r])):
            ecco_tile_number = ordered_ecco_tiles[r][c]
            face, min_row, min_col = ecco_tile_face_row_col_bounds(ecco_tile_number, llc, ecco_sNx,
                                                                   ecco_sNy)

            # print('  - Reading tile '+str(ecco_tile_number)+' from face '+str(face)+' (rows '+str(min_row)+' to '+str(min_row+ecco_sNy) + \
            #       ', cols '+str(min_col)+', '+str(min_col+ecco_sNx)+')')

            # print(ecco_tile_number,face,min_row,min_col)
            Runoff = runoff_faces[face][:, min_row:min_row + ecco_sNy, min_col:min_col + ecco_sNx]

            # rotate things as necessary
            for n in range(ordered_ecco_tile_rotations[r][c]):
                Runoff = np.rot90(Runoff,axes=(1,2))

            # put it into the big grid
            ecco_Runoff[:, r * ecco_sNy:(r + 1) * ecco_sNy, c * ecco_sNx:(c + 1) * ecco_sNx] = Runoff

    return(ecco_Runoff)

def read_ecco_exf_file(ecco_dir, file_prefix, year, llc):

    print('      - Reading from '+file_prefix+'_'+str(year))

    lon = np.arange(0,360,0.7031250)
    lon[lon>180] -= 360
    first_index = np.where(lon<0)[0][0]
    lon = np.concatenate([lon[first_index:],lon[:first_index]])

    del_lat = np.concatenate([np.array([0]),
                              np.array([0.6958694]),
                              np.array([0.6999817]),
                              np.array([0.7009048]),
                              np.array([0.7012634]),
                              np.array([0.7014313]),
                              0.7017418*np.ones((245,)),
                              np.array([0.7014313]),
                              np.array([0.7012634]),
                              np.array([0.7009048]),
                              np.array([0.6999817]),
                              np.array([0.6958694]),
                              ])
    lat = np.cumsum(del_lat)+-89.4628220

    grid = np.fromfile(os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'era_xx', file_prefix + '_' + str(year)), '>f4')
    n_timesteps = int(np.size(grid)/(len(lon)*len(lat)))
    grid = np.reshape(grid,(n_timesteps,len(lat),len(lon)))


    grid = np.concatenate([grid[:,:,first_index:],grid[:,:,:first_index]],axis=2)

    return(lon,lat,grid)


#################################################################################################
# pickup functions

def read_pickup_file_to_compact(pickup_file_path):

    Nr = 50
    print('      Reading from '+pickup_file_path)
    global_data, _, global_metadata = mds.rdmds(pickup_file_path, returnmeta=True)

    has_Nr = {'uvel': True, 'vvel': True, 'theta': True,
              'salt': True, 'gunm1': True, 'gvnm1': True,
              'gunm2': True, 'gvnm2': True, 'etan': False,
              'detahdt': False, 'etah': False}

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if has_Nr[var_name.strip().lower()]:
            end_row = start_row + Nr
        else:
            end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)

def read_tile_from_ECCO_compact(var_grid,tile_number, print_messages=False):

    sNx = 90
    sNy = 90

    if print_messages:
        print('Reading grid for tile number '+str(tile_number))

    # # adjust tile number to account for blank cells
    # tile_add = 6
    # if tile_number>2:
    #     tile_add+=1
    # if tile_number>4:
    #     tile_add+=1
    # if tile_number>23:
    #     tile_add+=4
    # if tile_number>26:
    #     tile_add+=2
    # tile_number += tile_add
    #
    # if print_messages:
    #     print('    - The tile number including blank cells is '+str(tile_number))

    # get the face number
    if tile_number < 28:
        face_number = 1
    if tile_number >= 28 and tile_number < 55:
        face_number = 2
    if tile_number >= 55 and tile_number<64:
        face_number = 3
    if tile_number >= 64 and tile_number<91:
        face_number = 4
    if tile_number >= 91:
        face_number = 5

    if print_messages:
        print('    - It is found in face '+str(face_number))

    face_indices_in_compact = {1:[0,9*sNx],
                               2:[9*sNx,(18)*sNx],
                               3:[(18)*sNx,(21)*sNx],
                               4:[(21)*sNx,(30)*sNx],
                               5:[(30)*sNx,(39)*sNx]}
    face_dimensions = {1:[9*sNy,3*sNx],
                       2:[9*sNy,3*sNx],
                       3:[3*sNy,3*sNx],
                       4:[3*sNy,9*sNx],
                       5:[3*sNy,9*sNx]}
    face_tile_dimensions = {1: [9, 3],
                           2: [9, 3],
                           3: [3, 3],
                           4: [3, 9],
                           5: [3, 9]}
    face_first_tile = {1:1,2:28,3:55,4:64,5:91}

    indices = face_indices_in_compact[face_number]
    dimensions = face_dimensions[face_number]

    if print_messages:
        print('    - The face will be dimension '+str(dimensions))
        print('    - The face will be read from row indices ' + str(indices)+' from the compact file')

    face_subset = var_grid[:,indices[0]:indices[1],:]
    face_subset = np.reshape(face_subset,(np.shape(face_subset)[0],dimensions[0],dimensions[1]))

    if print_messages:
        print('    - The tile number in the face is '+str((tile_number-face_first_tile[face_number]+1)))
        print('    - The number of tiles before the tile in the face are  '+str(face_tile_dimensions[face_number][1]))

    ll_row = sNy*int(((tile_number-face_first_tile[face_number]) // face_tile_dimensions[face_number][1]))
    ll_col = sNx*((tile_number-face_first_tile[face_number]) % face_tile_dimensions[face_number][1])

    if print_messages:
        print('    - The lower left row in the face is '+str(ll_row))
        print('    - The lower left col in the face is ' + str(ll_col))

    tile_subset = face_subset[:,ll_row:ll_row+sNy,ll_col:ll_col+sNx]

    # plt.imshow(tile_subset[0,:,:],origin='lower')
    # plt.title(str(tile_number-tile_add))
    # plt.show()

    return(tile_subset)

def read_ecco_pickup_to_stiched_grid(pickup_file_path,ordered_ecco_tiles,ordered_ecco_tile_rotations):

    ecco_Nr = 50
    ecco_sNx = 90
    ecco_sNy = 90

    # pickup_file_path = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init',pickup_file)
    var_names,row_bounds,compact_var_grids,global_metadata = read_pickup_file_to_compact(pickup_file_path)

    var_grids = []

    for vn in range(len(var_names)):
        compact_var_grid = compact_var_grids[vn]
        if var_names[vn].lower() not in ['etan', 'detahdt', 'etah']:
            ecco_grid = np.zeros((ecco_Nr, ecco_sNy * len(ordered_ecco_tiles), ecco_sNx * len(ordered_ecco_tiles[0])))
        else:
            ecco_grid = np.zeros((1, ecco_sNy * len(ordered_ecco_tiles), ecco_sNx * len(ordered_ecco_tiles[0])))

        for r in range(len(ordered_ecco_tiles)):
            ecco_tile_row = ordered_ecco_tiles[r]
            ecco_rotation_row = ordered_ecco_tile_rotations[r]
            for c in range(len(ordered_ecco_tiles[r])):

                # get the variable grid
                var_grid = read_tile_from_ECCO_compact(compact_var_grid, tile_number=ecco_tile_row[c])

                # rotate things as necessary
                for n in range(ecco_rotation_row[c]):
                    var_grid = np.rot90(var_grid, axes=(1, 2))

                # put it into the big grid
                ecco_grid[:, r * ecco_sNy:(r + 1) * ecco_sNy, c * ecco_sNx:(c + 1) * ecco_sNx] = var_grid

        var_grids.append(ecco_grid)

    return(var_names, var_grids, global_metadata)

def rotate_ecco_grids_to_natural_grids(var_names, var_grids, ecco_AngleCS, ecco_AngleSN):

    def rotate_velocity_vectors_to_natural(angle_cos, angle_sin, uvel, vvel):
        zonal_velocity = np.zeros_like(uvel)
        meridional_velocity = np.zeros_like(vvel)
        for k in range(np.shape(uvel)[0]):
            zonal_velocity[k,:,:] = angle_cos * uvel[k,:,:] - angle_sin * vvel[k,:,:]
            meridional_velocity[k,:,:] = angle_sin * uvel[k,:,:] + angle_cos * vvel[k,:,:]
        return (zonal_velocity, meridional_velocity)

    uvel_grid = var_grids[var_names.index('Uvel')]
    vvel_grid = var_grids[var_names.index('Vvel')]
    natural_uvel_grid, natural_vvel_grid = rotate_velocity_vectors_to_natural(ecco_AngleCS, ecco_AngleSN, uvel_grid, vvel_grid)
    var_grids[var_names.index('Uvel')] = natural_uvel_grid
    var_grids[var_names.index('Vvel')] = natural_vvel_grid

    gunm1_grid = var_grids[var_names.index('GuNm1')]
    gvnm1_grid = var_grids[var_names.index('GvNm1')]
    natural_gunm1_grid, natural_gvnm1_grid = rotate_velocity_vectors_to_natural(ecco_AngleCS, ecco_AngleSN, gunm1_grid, gvnm1_grid)
    var_grids[var_names.index('GuNm1')] = natural_gunm1_grid
    var_grids[var_names.index('GvNm1')] = natural_gvnm1_grid

    gunm2_grid = var_grids[var_names.index('GuNm2')]
    gvnm2_grid = var_grids[var_names.index('GvNm2')]
    natural_gunm2_grid, natural_gvnm2_grid = rotate_velocity_vectors_to_natural(ecco_AngleCS, ecco_AngleSN, gunm2_grid, gvnm2_grid)
    var_grids[var_names.index('GuNm2')] = natural_gunm2_grid
    var_grids[var_names.index('GvNm2')] = natural_gvnm2_grid

    return(var_grids)


#################################################################################################
# seaice pickup functions

def read_seaice_pickup_file_to_compact(pickup_file_path):

    print('      Reading from '+pickup_file_path)
    global_data, _, global_metadata = mds.rdmds(pickup_file_path, returnmeta=True)

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)

def read_ecco_seaice_pickup_to_stiched_grid(pickup_file_path,ordered_ecco_tiles,ordered_ecco_tile_rotations):

    ecco_sNx = 90
    ecco_sNy = 90

    # pickup_file_path = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init',pickup_file)
    var_names,row_bounds,compact_var_grids,global_metadata = read_seaice_pickup_file_to_compact(pickup_file_path)

    var_grids = []

    for vn in range(len(var_names)):
        compact_var_grid = compact_var_grids[vn]
        ecco_grid = np.zeros((1, ecco_sNy * len(ordered_ecco_tiles), ecco_sNx * len(ordered_ecco_tiles[0])))

        for r in range(len(ordered_ecco_tiles)):
            ecco_tile_row = ordered_ecco_tiles[r]
            ecco_rotation_row = ordered_ecco_tile_rotations[r]
            for c in range(len(ordered_ecco_tiles[r])):

                # get the variable grid
                var_grid = read_tile_from_ECCO_compact(compact_var_grid, tile_number=ecco_tile_row[c])

                # rotate things as necessary
                for n in range(ecco_rotation_row[c]):
                    var_grid = np.rot90(var_grid, axes=(1, 2))

                # put it into the big grid
                ecco_grid[:, r * ecco_sNy:(r + 1) * ecco_sNy, c * ecco_sNx:(c + 1) * ecco_sNx] = var_grid

        var_grids.append(ecco_grid)

    return(var_names, var_grids, global_metadata)

def rotate_ecco_seaice_grids_to_natural_grids(var_names, var_grids, ecco_AngleCS, ecco_AngleSN):

    def rotate_velocity_vectors_to_natural(angle_cos, angle_sin, uvel, vvel):
        zonal_velocity = np.zeros_like(uvel)
        meridional_velocity = np.zeros_like(vvel)
        for k in range(np.shape(uvel)[0]):
            zonal_velocity[k,:,:] = angle_cos * uvel[k,:,:] - angle_sin * vvel[k,:,:]
            meridional_velocity[k,:,:] = angle_sin * uvel[k,:,:] + angle_cos * vvel[k,:,:]
        return (zonal_velocity, meridional_velocity)

    uvel_grid = var_grids[var_names.index('siUICE')]
    vvel_grid = var_grids[var_names.index('siVICE')]
    natural_uvel_grid, natural_vvel_grid = rotate_velocity_vectors_to_natural(ecco_AngleCS, ecco_AngleSN, uvel_grid, vvel_grid)
    var_grids[var_names.index('siUICE')] = natural_uvel_grid
    var_grids[var_names.index('siVICE')] = natural_vvel_grid

    return(var_grids)

#################################################################################################
# ptracer pickup functions

def read_ptracer_pickup_file_to_compact(pickup_file_path):

    Nr = 50
    print('      Reading from '+pickup_file_path)
    global_data, _, global_metadata = mds.rdmds(pickup_file_path, returnmeta=True)

    var_names = []
    var_grids = []

    for vn in range(len(global_metadata['fldlist'])):
        var_grid = global_data[vn,:,:,:]
        var_grids.append(var_grid)
        var_names.append(global_metadata['fldlist'][vn].strip())

    return(var_names,var_grids,global_metadata)

def read_ecco_ptracer_pickup_to_stiched_grid(pickup_file_path,ordered_ecco_tiles,ordered_ecco_tile_rotations):

    ecco_Nr = 50
    ecco_sNx = 90
    ecco_sNy = 90

    # pickup_file_path = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init',pickup_file)
    var_names,compact_var_grids,global_metadata = read_ptracer_pickup_file_to_compact(pickup_file_path)

    var_grids = []

    for vn in range(len(var_names)):
        compact_var_grid = compact_var_grids[vn]
        ecco_grid = np.zeros((ecco_Nr, ecco_sNy * len(ordered_ecco_tiles), ecco_sNx * len(ordered_ecco_tiles[0])))

        for r in range(len(ordered_ecco_tiles)):
            ecco_tile_row = ordered_ecco_tiles[r]
            ecco_rotation_row = ordered_ecco_tile_rotations[r]
            for c in range(len(ordered_ecco_tiles[r])):

                # get the variable grid
                var_grid = read_tile_from_ECCO_compact(compact_var_grid, tile_number=ecco_tile_row[c])

                # rotate things as necessary
                for n in range(ecco_rotation_row[c]):
                    var_grid = np.rot90(var_grid, axes=(1, 2))

                # put it into the big grid
                ecco_grid[:, r * ecco_sNy:(r + 1) * ecco_sNy, c * ecco_sNx:(c + 1) * ecco_sNx] = var_grid

        var_grids.append(ecco_grid)

    return(var_names, var_grids, global_metadata)

#################################################################################################
# darwin pickup function

def read_darwin_pickup_file_to_compact(pickup_file_path):

    Nr = 50
    print('      Reading from '+pickup_file_path)
    global_data, _, global_metadata = mds.rdmds(pickup_file_path, returnmeta=True)

    # there is only one field in the darwin pickup

    var_names = [global_metadata['fldlist'][0]]
    var_grids = [global_data]

    return(var_names,var_grids,global_metadata)

def read_ecco_darwin_pickup_to_stiched_grid(pickup_file_path,ordered_ecco_tiles,ordered_ecco_tile_rotations):

    ecco_Nr = 50
    ecco_sNx = 90
    ecco_sNy = 90

    # pickup_file_path = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init',pickup_file)
    var_names,compact_var_grids,global_metadata = read_darwin_pickup_file_to_compact(pickup_file_path)

    var_grids = []

    for vn in range(len(var_names)):
        compact_var_grid = compact_var_grids[vn]
        ecco_grid = np.zeros((ecco_Nr, ecco_sNy * len(ordered_ecco_tiles), ecco_sNx * len(ordered_ecco_tiles[0])))

        for r in range(len(ordered_ecco_tiles)):
            ecco_tile_row = ordered_ecco_tiles[r]
            ecco_rotation_row = ordered_ecco_tile_rotations[r]
            for c in range(len(ordered_ecco_tiles[r])):

                # get the variable grid
                var_grid = read_tile_from_ECCO_compact(compact_var_grid, tile_number=ecco_tile_row[c])

                # rotate things as necessary
                for n in range(ecco_rotation_row[c]):
                    var_grid = np.rot90(var_grid, axes=(1, 2))

                # put it into the big grid
                ecco_grid[:, r * ecco_sNy:(r + 1) * ecco_sNy, c * ecco_sNx:(c + 1) * ecco_sNx] = var_grid

        var_grids.append(ecco_grid)

    return(var_names, var_grids, global_metadata)

#################################################################################################
# VELMASS to VEL conversion

def convert_velmass_to_vel(ecco_dir, llc, var_name, year):
    # to convert to UVEL and VVEL, we will use the relationship explained here:
    # https://ecco-v4-python-tutorial.readthedocs.io/ECCO_v4_Volume_budget_closure.html

    #######################################################################################
    # read in the fields

    # get the velmass file
    ds = nc4.Dataset(os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','nctiles_monthly', var_name + 'MASS', var_name + 'MASS_' + str(year) + '.nc'))
    if var_name == 'UVEL':
        i = ds.variables['i_g'][:]
        j = ds.variables['j'][:]
    if var_name == 'VVEL':
        i = ds.variables['i'][:]
        j = ds.variables['j_g'][:]
    k = ds.variables['k'][:]
    tile_var = ds.variables['tile'][:]
    time = ds.variables['time'][:]
    timestep_var = ds.variables['timestep'][:]
    velmass = ds.variables[var_name + 'MASS'][:, :, :, :, :]
    ds.close()

    # get the ETAN field
    ds = nc4.Dataset(os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','nctiles_monthly','ETAN','ETAN_' + str(year) + '.nc'))
    etan = ds.variables['ETAN'][:, :, :, :]
    ds.close()

    if var_name=='UVEL':
        file_path = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init','hFacW.data')
    if var_name=='VVEL':
        file_path = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init','hFacS.data')
    hFac_faces = read_ecco_field_to_faces(file_path, llc, dim=3)
    hFac_tiles = read_ecco_faces_to_tiles(hFac_faces, llc, dim=3)

    depth_file_path = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'input_init', 'bathy270_filled_noCaspian_r4')
    depth_faces = read_ecco_field_to_faces(depth_file_path, llc, dim=2)
    depth_tiles = read_ecco_faces_to_tiles(depth_faces, llc, dim=2)

    #######################################################################################
    # do the conversion

    vel = np.zeros_like(velmass)
    for timestep in range(np.shape(velmass)[0]):
        for tile in range(np.shape(velmass)[2]):
            s_star = 1 + etan[timestep, tile, :, :] / depth_tiles[tile+1][:, :]
            vel[timestep, :, tile, :, :] = velmass[timestep, :, tile, :, :] / (hFac_tiles[tile+1] * s_star)

    # double check its masked
    for timestep in range(np.shape(velmass)[0]):
        for tile in range(np.shape(velmass)[2]):
            subset = vel[timestep, :, tile, :, :]
            subset[hFac_tiles[tile+1][:, :, :] == 0] = 0
            vel[timestep, :, tile, :, :] = subset

    # #######################################################################################
    # # maybe plot a sanity check
    # plt.subplot(1,2,1)
    # C = plt.imshow(velmass[3,0,7,:,:],origin='lower')
    # plt.colorbar(C)
    # plt.title('velmass')
    # plt.subplot(1, 2, 2)
    # C = plt.imshow(vel[3, 0, 7, :, :], origin='lower')
    # plt.colorbar(C)
    # plt.title('vel')
    # plt.show()

    #######################################################################################
    # write out the field
    if var_name not in os.listdir(os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','nctiles_monthly')):
        os.mkdir(os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','nctiles_monthly',var_name))

    ds = nc4.Dataset(os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','nctiles_monthly', var_name, var_name + '_' + str(year) + '.nc'), 'w')

    if var_name=='UVEL':
        ds.createDimension('i_g', np.size(i))
        ds.createDimension('j', np.size(j))
    if var_name=='VVEL':
        ds.createDimension('i', np.size(i))
        ds.createDimension('j_g', np.size(j))
    ds.createDimension('k', np.size(k))
    ds.createDimension('tile', np.size(tile_var))
    ds.createDimension('time', np.size(time))

    if var_name=='UVEL':
        var = ds.createVariable('i_g', 'f4', ('i_g',))
        var[:] = i
        var = ds.createVariable('j', 'f4', ('j',))
        var[:] = j
    if var_name=='VVEL':
        var = ds.createVariable('i', 'f4', ('i',))
        var[:] = i
        var = ds.createVariable('j_g', 'f4', ('j_g',))
        var[:] = j

    var = ds.createVariable('k', 'f4', ('k',))
    var[:] = k

    var = ds.createVariable('tile', 'f4', ('tile',))
    var[:] = tile_var

    var = ds.createVariable('time', 'f4', ('time',))
    var[:] = time

    var = ds.createVariable('timestep', 'f4', ('time',))
    var[:] = timestep_var

    if var_name=='UVEL':
        var = ds.createVariable('UVEL', 'f4', ('time', 'k', 'tile', 'i_g', 'j'))
        var[:, :, :, :, :] = vel
    if var_name=='VVEL':
        var = ds.createVariable('VVEL', 'f4', ('time', 'k', 'tile', 'i', 'j_g'))
        var[:, :, :, :, :] = vel

    ds.close()





