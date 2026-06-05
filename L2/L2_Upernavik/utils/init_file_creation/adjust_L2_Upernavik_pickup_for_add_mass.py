
import os
import numpy as np
import matplotlib.pyplot as plt
from MITgcmutils import mds
import netCDF4 as nc4
from datetime import timedelta, datetime
import argparse

def read_grid_geometry_from_nc(config_dir,model_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    hFacC = ds.variables['HFacC'][:, :, :]
    delR = np.array(ds.variables['drF'][:])
    ds.close()
    return(XC,YC,hFacC,delR)

def read_pickup_file_to_compact(pickup_file_path, Nr):

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

def iter_number_to_date(iter_number,model_level):
    a=1
    if model_level=='L0':
        seconds_per_iter = 1200
    if model_level=='L1':
        seconds_per_iter = 150
    if model_level=='L2':
        seconds_per_iter = 60

    total_seconds = iter_number*seconds_per_iter
    date = datetime(1992,1,1) + timedelta(seconds=total_seconds)
    return(date)

def stack_grids_to_pickup(interp_grids):
    counter = 0
    for grid in interp_grids:

        if counter == 0:
            pickup_grid = grid
        else:
            pickup_grid = np.concatenate([pickup_grid, grid], axis=0)

        counter += 1
    return(pickup_grid)

def write_pickup_file(output_file,dtype,pickup_grid,subset_metadata):

    # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+"\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [   "+str(subset_metadata['nrecords'][0])+" ];\n"
    output += " timeStepNumber = [ "+"{:10d}".format(subset_metadata['timestepnumber'][0])+" ];\n"
    time_interval_exponent = int(np.log10(subset_metadata['timeinterval'][0][0]))
    time_interval_base = subset_metadata['timeinterval'][0][0] / (10 ** time_interval_exponent)
    output += " timeInterval = [  "+"{:.12f}".format(time_interval_base) + "E+" + "{:02d}".format(time_interval_exponent)+  " ];\n"
    output += " nFlds = [   "+str(subset_metadata['nflds'][0])+" ];\n"
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

def adjust_pickup(config_dir, model_name, pickup_iteration):
    XC, YC, hFacC, delR = read_grid_geometry_from_nc(config_dir,model_name)
    Nr = len(delR)

    date = iter_number_to_date(pickup_iteration, model_level='L2')
    year = date.year

    add_mass = True

    pickup_file_path = os.path.join(config_dir,'L2',model_name,'run_baseline','pickup.'+'{:010d}'.format(pickup_iteration))

    var_names, row_bounds, var_grids, global_metadata = read_pickup_file_to_compact(pickup_file_path, Nr)

    masked_var_grids = []
    for v in range(len(var_names)):
        if np.shape(var_grids[v])[0] == 1:
            new_var_grid = var_grids[v]# * mask[0, :, :]
        else:
            new_var_grid = var_grids[v]# * mask
        masked_var_grids.append(new_var_grid)

        if var_names[v]=='GvNm2' and add_mass and 'addMass' not in var_names:
            print('    - Adding a zero add mass field to the pickup')
            var_names.insert(var_names.index('GvNm2')+1, 'AddMass')
            masked_var_grids.append(np.zeros_like(new_var_grid))

            global_metadata['fldlist'] = var_names
            global_metadata['nflds'][0] = global_metadata['nflds'][0]+1

    pickup_grid = stack_grids_to_pickup(masked_var_grids)
    print(np.shape(pickup_grid))
    global_metadata['nrecords'][0] = np.shape(pickup_grid)[0]

    output_file = os.path.join(config_dir, 'L2', model_name, 'run_baseline_iceplume',
                               'pickup.' + '{:010d}'.format(pickup_iteration)+'_addMass')
    write_pickup_file(output_file, '>f8', pickup_grid, global_metadata)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)
    args = parser.parse_args()

    config_dir = args.config_dir

    model_name = 'L2_Upernavik'

    pickup_iteration = 30332160

    adjust_pickup(config_dir, model_name, pickup_iteration)












