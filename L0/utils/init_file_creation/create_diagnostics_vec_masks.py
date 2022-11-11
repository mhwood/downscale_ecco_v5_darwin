import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
from math import radians, cos, sin, asin, sqrt
import simplegrid as sg
from scipy.interpolate import griddata
from ecco_v4_py.llc_array_conversion import llc_faces_to_compact, llc_compact_to_faces
import argparse
import ast


def read_global_XC_YC(ecco_dir,llc):

    grid_file_dir = os.path.join(ecco_dir,'mitgrid_tiles')
    # grid_file_dir = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','mitgrid_tiles')
    XC_faces = {}
    YC_faces = {}
    for i in range(1,6):
        if i<3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc, 3*llc)
        if i==3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc, llc)
        if i>3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), 3*llc, llc)
        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face

    input_init_dir = os.path.join(ecco_dir, 'input_init')
    bathy_compact = np.fromfile(os.path.join(input_init_dir, 'bathy_llc' + str(llc)), '>f4')
    bathy_compact = np.reshape(bathy_compact, (13 * llc, llc))
    bathy_faces = llc_compact_to_faces(bathy_compact, less_output=True)

    return(XC_faces,YC_faces, bathy_faces)

def read_model_geometry_from_nc(config_dir,model_name):
    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',model_name+'_grid.nc'))
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    ds.close()
    return(XC, YC)

def great_circle_distance(lon_ref, lat_ref, Lon, Lat):
    earth_radius = 6371000
    lon_ref_radians = np.radians(lon_ref)
    lat_ref_radians = np.radians(lat_ref)
    lons_radians = np.radians(Lon)
    lats_radians = np.radians(Lat)
    lat_diff = lats_radians - lat_ref_radians
    lon_diff = lons_radians - lon_ref_radians
    d = np.sin(lat_diff * 0.5) ** 2 + np.cos(lat_ref_radians) * np.cos(lats_radians) * np.sin(lon_diff * 0.5) ** 2
    h = 2 * earth_radius * np.arcsin(np.sqrt(d))
    return(h)

def subset_geometry_to_boundary(boundary, XC, YC):

    # subset to the boundary
    if boundary=='south':
        XC_subset = XC[:1, :]
        YC_subset = YC[:1, :]

    if boundary=='west':
        XC_subset = XC[:, :1]
        YC_subset = YC[:, :1]

    if boundary=='north':
        XC_subset = XC[-1:, :]
        YC_subset = YC[-1:, :]

    if boundary=='east':
        XC_subset = XC[:, -1:]
        YC_subset = YC[:, -1:]

    return(XC_subset, YC_subset)

def create_boundary_masks(llc, resolution, XC_faces, YC_faces, bathy_faces,
                         northern_points, southern_points, eastern_points, western_points):

    all_masks = []
    mask_dicts = []
    n_mask_points = []
    boundaries = ['north','south','east','west']

    for b in range(len(boundaries)):
        print('        - Generating the '+boundaries[b]+' mask')
        if b==0:
            points = northern_points
        if b==1:
            points = southern_points
        if b ==2:
            points = eastern_points
        if b==3:
            points = western_points


        if len(points)>0:
            mask_faces = {}
            mask_faces[1] = np.zeros((3 * llc, llc))
            mask_faces[2] = np.zeros((3 * llc, llc))
            mask_faces[3] = np.zeros((llc, llc))
            mask_faces[4] = np.zeros((llc,3 * llc))
            mask_faces[5] = np.zeros((llc,3 * llc))

            mask_dict = np.zeros((llc * llc * 5, 3))

            counter = 1
            for face in [1,2,3,4,5]:
                XC = XC_faces[face]
                YC = YC_faces[face]
                bathy = bathy_faces[face]
                for n in range(np.shape(points)[0]):
                    x = points[n,0]
                    y = points[n,1]
                    # if n == int(0.1*np.shape(points)[0]):
                    #     print('            - face '+str(face)+' is 10% complete')
                    if n == int(0.2*np.shape(points)[0]):
                        print('            - face '+str(face)+' is 20% complete')
                    # if n == int(0.3*np.shape(points)[0]):
                    #     print('            - face '+str(face)+' is 30% complete')
                    if n == int(0.4*np.shape(points)[0]):
                        print('            - face '+str(face)+' is 40% complete')
                    # if n == int(0.5*np.shape(points)[0]):
                    #     print('            - face '+str(face)+' is 50% complete')
                    if n == int(0.6*np.shape(points)[0]):
                        print('            - face '+str(face)+' is 60% complete')
                    # if n == int(0.7*np.shape(points)[0]):
                    #     print('            - face '+str(face)+' is 70% complete')
                    if n == int(0.8*np.shape(points)[0]):
                        print('            - face '+str(face)+' is 80% complete')
                    # if n == int(0.9*np.shape(points)[0]):
                    #     print('            - face '+str(face)+' is 90% complete')
                    if boundaries[b]=='surface':
                        dist = ((XC-x)**2 + (YC-y)**2)**0.5
                    else:
                        dist = great_circle_distance(x, y, XC, YC)
                    # print('using the new formula')
                    rows, cols = np.where(dist <= resolution*np.sqrt(2))
                    if len(rows)>0:
                        for i in range(len(rows)):
                            row = rows[i]
                            col = cols[i]
                            if bathy[row,col]<0:
                                if mask_faces[face][row,col]==0:
                                    mask_faces[face][row,col] = counter
                                    mask_dict[counter-1,0] = face
                                    mask_dict[counter-1,1] = row
                                    mask_dict[counter-1,2] = col
                                    counter += 1
            # print(counter-1)
            n_points = np.sum(np.sum(mask_dict[:, 0] != 0))
            n_mask_points.append(n_points)
            if n_points>0:
                mask_dict = mask_dict[:np.sum(mask_dict[:, 0] != 0), :]
                print('            - This mask will have '+str(counter)+' points')
                all_masks.append(mask_faces)
                mask_dicts.append(mask_dict)
            else:
                all_masks.append({})
                mask_dicts.append([])
        else:
            all_masks.append({})
            mask_dicts.append([])

    return(all_masks, mask_dicts, n_mask_points, boundaries)

def create_surface_mask(llc, XC_faces, YC_faces, bathy_faces, ecco_sNx, ecco_sNy,
                       ecco_surface_tile_faces, ecco_surface_tile_min_rows, ecco_surface_tile_min_cols):

    mask_faces = {}
    mask_faces[1] = np.zeros((3 * llc, llc))
    mask_faces[2] = np.zeros((3 * llc, llc))
    mask_faces[3] = np.zeros((llc, llc))
    mask_faces[4] = np.zeros((llc,3 * llc))
    mask_faces[5] = np.zeros((llc,3 * llc))

    point_buffer = 2

    mask_dict = np.zeros((ecco_sNx * ecco_sNy * len(ecco_surface_tile_faces), 3))

    counter = 1
    for f in range(len(ecco_surface_tile_faces)):
        face = ecco_surface_tile_faces[f]
        XC = XC_faces[face]
        YC = YC_faces[face]
        bathy = bathy_faces[face]

        for row in range(ecco_surface_tile_min_rows[f] - point_buffer,ecco_surface_tile_min_rows[f]+ecco_sNy):
            for col in range(ecco_surface_tile_min_cols[f], ecco_surface_tile_min_cols[f] + ecco_sNx + point_buffer):
                if bathy[row,col]<0:
                    if mask_faces[face][row,col]==0:
                        mask_faces[face][row, col] = counter
                        mask_dict[counter - 1, 0] = face
                        mask_dict[counter - 1, 1] = row
                        mask_dict[counter - 1, 2] = col
                        counter += 1

    # print(counter-1)
    mask_dict = mask_dict[:np.sum(mask_dict[:, 0] != 0), :]

    print('            - This mask will have '+str(counter)+' points')

    return(mask_faces, mask_dict)


def output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names):
    if output_file_name in os.listdir(output_dir):
        os.remove(os.path.join(output_dir,output_file_name))

    ds = nc4.Dataset(os.path.join(output_dir,output_file_name),'w')

    for m in range(len(mask_names)):
        print(mask_names[m],np.shape(all_mask_dicts[m]))
        if len(all_mask_dicts[m])>0:
            grp = ds.createGroup(mask_names[m])
            grp.createDimension('n_points', np.shape(all_mask_dicts[m])[0])
            var = grp.createVariable('source_faces', 'i4', ('n_points',))
            var[:] = all_mask_dicts[m][:, 0].astype(int)
            var = grp.createVariable('source_rows', 'i4', ('n_points',))
            var[:] = all_mask_dicts[m][:, 1].astype(int)
            var = grp.createVariable('source_cols', 'i4', ('n_points',))
            var[:] = all_mask_dicts[m][:, 2].astype(int)

    ds.close()



########################################################################################################################

def create_dv_masks(config_dir, ecco_path,L1_model_names,print_status):

    if print_status:
        print('Creating the diagnostics_vec masks to use in the L0 domain')

    llc = 270
    resolution = 25e3 # in m

    if print_status:
        print('    - Reading in the L0 domain files')

    # read the mitgrids to faces
    ecco_dir = os.path.join(ecco_path,'LLC'+str(llc)+'_Files')
    XC_faces, YC_faces, bathy_faces = read_global_XC_YC(ecco_dir, llc)

    all_mask_names = []
    all_masks = []
    all_mask_dicts = []

    for L1_model_name in L1_model_names:

        ###############################################################################################
        # Read in the L1 model information

        XC, YC = read_model_geometry_from_nc(config_dir,L1_model_name)

        southern_points = np.column_stack([XC[:1, :].ravel(), YC[:1, :].ravel()])
        northern_points = np.column_stack([XC[-1:, :].ravel(), YC[-1:, :].ravel()])
        western_points = np.column_stack([XC[:, :1].ravel(), YC[:, :1].ravel()])
        eastern_points = np.column_stack([XC[:, -1:].ravel(), YC[:, -1:].ravel()])

        # ###############################################################################################
        # # Mask masks of the boundary points in the faces
        #
        # print('    - Generating the surface mask')
        # surface_mask, surface_mask_dict = create_surface_mask(llc, XC_faces, YC_faces, bathy_faces, ecco_sNx, ecco_sNy,
        #                                                       ecco_surface_tile_faces, ecco_surface_tile_min_rows,
        #                                                       ecco_surface_tile_min_cols)
        # all_mask_names += [L1_model_name + '_surface']
        # all_masks += [surface_mask]
        # all_mask_dicts += [surface_mask_dict]

        ###############################################################################################
        # Mask masks of the boundary points in the faces

        print('    - Generating the BC masks')
        all_mask_names += [L1_model_name+'_north',L1_model_name+'_south',L1_model_name+'_east',L1_model_name+'_west']
        BC_masks, BC_mask_dicts, n_mask_points, boundaries = create_boundary_masks(llc, resolution, XC_faces, YC_faces, bathy_faces,
                                       northern_points, southern_points, eastern_points, western_points)
        all_masks += BC_masks
        all_mask_dicts += BC_mask_dicts

    ###############################################################################################
    # Convert the masks to compact and output as a binary file

    for m in range(len(all_mask_names)):
        mask_faces = all_masks[m]
        if bool(set(list(mask_faces.keys())) & set(list([1,2,3,4,5]))):
            print('    - Outputting the '+all_mask_names[m]+' mask to binary')
            compact_mask = llc_faces_to_compact(mask_faces)
            output_name = all_mask_names[m]
            if len('dv/'+output_name)>30:
                raise ValueError('diagnostics_vec cannot take names longer than 30 chars')
            output_file = os.path.join(config_dir, 'L0', 'input', all_mask_names[m] + '_mask.bin')
            compact_mask.ravel('C').astype('>f4').tofile(output_file)

    output_dir = os.path.join(config_dir, 'L0', 'input')
    output_file_name = 'L0_dv_mask_reference_dict.nc'
    output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,all_mask_names)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The name of the configuration.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=True)

    parser.add_argument("-m", "--model_names", action="store",
                        help="List of model names to create masks for (e.g. L1_W_Greenland, L1_GOM)", dest="model_names",
                        type=str, required=True, nargs='+')

    parser.add_argument("-p", "--print_status", action="store",
                        help="Print status of routine (1 for True, 0 for False).", dest="print_status",
                        type=int, required=False, default=1)

    args = parser.parse_args()
    ecco_path = args.ecco_path
    config_dir = args.config_dir
    L1_model_names = args.model_names
    print_status = args.print_status

    if print_status>0:
        print_status=True
    else:
        print_status=False

    create_dv_masks(config_dir,ecco_path,L1_model_names,print_status)
