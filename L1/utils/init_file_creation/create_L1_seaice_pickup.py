
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
from MITgcmutils import mds
from scipy.interpolate import griddata
import sys


def read_grid_geometry(config_dir,model_name):

    file_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]
    delR = ds.variables['drF'][:]
    ds.close()

    return(XC, YC, AngleCS, AngleSN, delR)

def read_grid_mask(config_dir,model_name,var_name):

    file_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(file_path)

    if var_name in ['siVICE']:
        hFac = 'S'
    elif var_name in ['siUICE']:
        hFac = 'W'
    else:
        hFac = 'C'

    HFac = ds.variables['HFac' + hFac][:, :, :]
    if hFac == 'W':
        HFac = HFac[:, :, :-1]
    if hFac == 'S':
        HFac = HFac[:, :-1, :]

    ds.close()

    mask = np.copy(HFac)
    mask[mask>0]=1

    return(mask)

def read_interpolation_grid(df,config_dir, model_name, var_name,
                            ecco_XC, ecco_YC, ecco_wet_cells, ecco_wet_cells_on_tile_domain, XC, YC, print_level):

    interpolation_grid_file = model_name+'_interpolation_grid.nc'

    if interpolation_grid_file not in os.listdir(os.path.join(config_dir,'nc_grids')):

        if print_level >= 3:
            print('            - Creating the interpolation grids for this model (hasnt been done yet!)')

        domain_wet_cells_3D_C = read_grid_mask(config_dir,model_name,'Theta')
        domain_wet_cells_3D_S = read_grid_mask(config_dir, model_name, 'Vvel')
        domain_wet_cells_3D_W = read_grid_mask(config_dir, model_name, 'Uvel')

        interpolation_type_grid_C, source_row_grid_C, source_col_grid_C, source_level_grid_C, \
        interpolation_type_grid_S, source_row_grid_S, source_col_grid_S, source_level_grid_S, \
        interpolation_type_grid_W, source_row_grid_W, source_col_grid_W, source_level_grid_W = \
            df.create_interpolation_grids(ecco_XC, ecco_YC, ecco_wet_cells, ecco_wet_cells_on_tile_domain,
                               XC, YC, domain_wet_cells_3D_C, domain_wet_cells_3D_S, domain_wet_cells_3D_W)

        df.write_interpolation_grid_to_nc(config_dir, model_name,
                                       interpolation_type_grid_C, source_row_grid_C, source_col_grid_C, source_level_grid_C,
                                       interpolation_type_grid_S, source_row_grid_S, source_col_grid_S, source_level_grid_S,
                                       interpolation_type_grid_W, source_row_grid_W, source_col_grid_W, source_level_grid_W)

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',interpolation_grid_file))
    if var_name in ['Vvel','SIvice','siVICE','VWIND','VSTRESS']:
        grp = ds.groups['S']
    elif var_name in ['Uvel','SIuice','siUICE','UWIND','USTRESS']:
        grp = ds.groups['W']
    else:
        grp = ds.groups['C']

    interpolation_type_grid = grp.variables['interp_type'][:,:,:]
    source_rows_grid = grp.variables['source_rows'][:, :, :]
    source_cols_grid = grp.variables['source_cols'][:, :, :]
    source_levels_grid = grp.variables['source_levels'][:,:,:]

    ds.close()

    return(interpolation_type_grid, source_rows_grid, source_cols_grid, source_levels_grid)

def interpolate_ecco_wetgrid_to_domain(XC, YC, ecco_XC, ecco_YC, ecco_wet_cells):

    ecco_wet_cells_on_tile_domain = np.zeros((1,np.shape(XC)[0],np.shape(XC)[1]))

    points = np.hstack([np.reshape(ecco_XC, (np.size(ecco_XC), 1)),
                        np.reshape(ecco_YC, (np.size(ecco_YC), 1))])

    values = np.reshape(ecco_wet_cells[0,:,:], (np.size(ecco_wet_cells[0,:,:]), 1))
    grid = griddata(points, values, (XC, YC), method='nearest')[:, :, 0]
    ecco_wet_cells_on_tile_domain[0,:,:] = grid

    return(ecco_wet_cells_on_tile_domain)

def rotate_interpolated_grids_to_domain(var_names, interp_grids, AngleCS, AngleSN):

    def rotate_velocity_vectors_to_domain(angle_cos, angle_sin, zonal_vel, meridional_vel):
        uvel = np.zeros_like(zonal_vel)
        vvel = np.zeros_like(meridional_vel)
        for k in range(np.shape(uvel)[0]):
            uvel[k, : ,:] = angle_cos * zonal_vel[k,:,:] + angle_sin * meridional_vel[k,:,:]
            vvel[k, : ,:] = -1 * angle_sin * zonal_vel[k,:,:] + angle_cos * meridional_vel[k,:,:]
        return (uvel, vvel)

    uvel_grid_index = var_names.index('siUICE')
    vvel_grid_index = var_names.index('siVICE')
    natural_uvel = interp_grids[uvel_grid_index]
    natural_vvel = interp_grids[vvel_grid_index]
    uvel, vvel = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                   natural_uvel, natural_vvel)
    interp_grids[uvel_grid_index] = uvel
    interp_grids[vvel_grid_index] = vvel

    return(interp_grids)

def stack_grids_to_pickup(interp_grids, var_names):

    for i in range(len(interp_grids)):
        print('       - Stacking field for '+var_names[i])
        var_grid = interp_grids[i]
        if i==0:
            pickup_grid = var_grid
        else:
            pickup_grid = np.concatenate([pickup_grid,var_grid], axis=0)

    return(pickup_grid)

def write_seaice_pickup_file(output_file,dtype,pickup_grid,subset_metadata):

    # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+"\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [    "+str(subset_metadata['nrecords'][0])+" ];\n"
    output += " timeStepNumber = [ "+"{:10d}".format(subset_metadata['timestepnumber'][0])+" ];\n"
    time_interval_exponent = int(np.log10(subset_metadata['timeinterval'][0][0]))
    time_interval_base = subset_metadata['timeinterval'][0][0] / (10 ** time_interval_exponent)
    output += " timeInterval = [  "+"{:.12f}".format(time_interval_base) + "E+" + "{:02d}".format(time_interval_exponent)+  " ];\n"
    output += " nFlds = [    "+str(subset_metadata['nflds'][0])+" ];\n"
    output += " fldList = {\n "
    for var_name in subset_metadata['fldlist']:
        output += "'"+var_name
        for i in range(8-len(var_name)):
            output+= " "
        output+="' "
    output += "\n };"

    f = open(output_file+'.meta','w')
    f.write(output)
    f.close()

def create_L1_seaice_pickup_file(config_dir, model_name,
                                      ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                                      parent_model_pickup_iteration, print_level, use_interpolation_grids = True):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import ecco_functions as ef

    if print_level >= 1:
        print('    - Creating the pickup file for the ' + model_name + ' model from ECCO data')

    if print_level >= 1:
        print('    - Reading in the L1 tile geometry')
    # step 0: get the model domain
    XC, YC, AngleCS, AngleSN, delR = read_grid_geometry(config_dir, model_name)
    Nr = len(delR)

    if print_level>=1:
        print('    - Reading in the ECCO tile geometry')
    # step 1: get the ecco faces geometry
    # ecco_AngleCS, ecco_AngleSN, ecco_hfacC,
    ecco_XC, ecco_YC, ecco_AngleCS, ecco_AngleSN, ecco_hfacC, _ = \
        ef.read_ecco_grid_geometry(ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations)

    if model_name == 'L1_mac_delta':
        ecco_XC += 360


    # plt.plot(ecco_XC, ecco_YC, 'b.')
    # plt.plot(XC, YC, 'k.')
    # plt.show()

    # C = plt.imshow(ecco_hfacC[0,:,:], origin='lower')
    # plt.colorbar(C)
    # plt.show()

    if print_level >= 1:
        print('    - Reading in the ECCO pickup file')
    pickup_file = 'pickup_seaice.'+'{:010d}'.format(parent_model_pickup_iteration)
    pickup_file_path = os.path.join(config_dir, 'L0', 'run', pickup_file)
    var_names, var_grids, global_metadata = ef.read_ecco_seaice_pickup_to_stiched_grid(pickup_file_path, ordered_ecco_tiles, ordered_ecco_tile_rotations)

    # C = plt.imshow(var_grids[4][0,:,:], origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
    # plt.colorbar(C)
    # plt.show()

    if print_level>=1:
        print('    - Rotating oriented fields to natural coordinates')
    var_grids = ef.rotate_ecco_seaice_grids_to_natural_grids(var_names, var_grids, ecco_AngleCS, ecco_AngleSN)

    # C = plt.imshow(var_grids[1][0,:,:], origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
    # plt.colorbar(C)
    # plt.show()

    ecco_wet_cells = np.copy(ecco_hfacC)
    ecco_wet_cells[ecco_wet_cells>0]=1

    # make some bins where the tiles will be stored
    interp_grids = []

    ####################################################################################################################


    ecco_wet_cells_on_tile_domain = interpolate_ecco_wetgrid_to_domain(XC, YC, ecco_XC, ecco_YC, ecco_wet_cells)

    for vn in range(len(var_names)):
        var_name = var_names[vn]

        if var_name not in []:  # used for testing
            ecco_grid = var_grids[var_names.index(var_name)]
            if print_level >= 3:
                print('            - Downscaling ' + var_name)

            domain_wet_cells_3D = read_grid_mask(config_dir, model_name, var_name)
            domain_wet_cells_3D = domain_wet_cells_3D[:1, :, :]

            if use_interpolation_grids:
                L1_interpolation_mask, L1_source_rows, L1_source_cols, L1_source_levels = \
                    read_interpolation_grid(df, config_dir, model_name, var_name,
                                            ecco_XC, ecco_YC, ecco_wet_cells, ecco_wet_cells_on_tile_domain,
                                            XC, YC, print_level)

            # plt.subplot(1, 3, 1)
            # plt.imshow(ecco_wet_cells[0, :, :], origin='lower')
            # plt.subplot(1, 3, 2)
            # plt.imshow(ecco_wet_cells_on_tile_domain[0, :, :], origin='lower')
            # plt.subplot(1, 3, 3)
            # plt.imshow(domain_wet_cells_3D[0, :, :], origin='lower')
            # plt.show()

            mean_vertical_difference = 0
            subset_copy = np.copy(ecco_grid)

            # plt.subplot(1,2,1)
            # plt.imshow(subset_copy[:,10,:])
            # plt.subplot(1, 2, 2)
            # plt.imshow(ecco_grid[:, 10, :])
            # plt.show()

            if print_level >= 4:
                printing = True
            else:
                printing = False

            if use_interpolation_grids:
                interp_field = df.downscale_3D_field_with_interpolation_mask(ecco_XC, ecco_YC,
                                                                             ecco_grid, ecco_wet_cells,
                                                                             ecco_wet_cells_on_tile_domain,
                                                                             XC, YC,
                                                                             domain_wet_cells_3D,
                                                                             L1_interpolation_mask, L1_source_rows,
                                                                             L1_source_cols, L1_source_levels,
                                                                             printing=printing)
            else:
                interp_field = df.downscale_3D_field_with_zeros(ecco_XC, ecco_YC,
                                                     ecco_grid, ecco_wet_cells,
                                                     ecco_wet_cells_on_tile_domain,
                                                     XC, YC, domain_wet_cells_3D,
                                                     mean_vertical_difference=0, fill_downward=True,
                                                     printing=False)

            if np.sum(np.isnan(interp_field))>0:
                print('Setting '+str(np.sum(np.isnan(interp_field)))+' values to 0 in this grid')
                interp_field[np.isnan(interp_field)] = 0

            # plt.subplot(2, 3, 1)
            # plt.imshow(ecco_wet_cells[0, :, :], origin='lower')
            # plt.subplot(2, 3, 2)
            # plt.imshow(ecco_wet_cells_on_tile_domain[0, :, :], origin='lower')
            # plt.subplot(2, 3, 3)
            # plt.imshow(domain_wet_cells_3D[0, :, :], origin='lower')
            # plt.subplot(2, 3, 4)
            # plt.imshow(ecco_grid[0, :, :], origin='lower')
            # plt.subplot(2, 3, 5)
            # plt.imshow(interp_field[0, :, :], origin='lower')
            # plt.show()

            interp_grids.append(interp_field)


    interp_grids = rotate_interpolated_grids_to_domain(var_names, interp_grids, AngleCS, AngleSN)

    # for face in faces:
    #     plt.imshow(interp_grid_faces[0][face][0,:,:],origin='lower')
    #     plt.title(var_names[0]+' - face '+str(face))
    #     plt.show()

    print('    - Stacking the interpolated fields into a compact pickup grid')
    pickup_grid = stack_grids_to_pickup(interp_grids,var_names)

    print('    - Outputting the compact pickup grid to the input directory')
    pickup_metadata = dict(global_metadata)
    output_dir = os.path.join(config_dir, 'L1', model_name, 'input')
    output_file = os.path.join(output_dir, 'pickup_seaice.' + '{:010d}'.format(4*parent_model_pickup_iteration))
    dtype = '>f8'
    pickup_metadata['timestepnumber'] = [4*int(pickup_metadata['timestepnumber'][0])]
    # pickup_metadata['nrecords'] = [np.shape(pickup_grid)[0]]
    # pickup_metadata['fldlist'] = var_names

    write_seaice_pickup_file(output_file, dtype, pickup_grid, pickup_metadata)
