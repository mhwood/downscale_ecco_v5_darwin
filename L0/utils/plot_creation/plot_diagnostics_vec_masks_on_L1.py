import os
import numpy as np
import simplegrid as sg
import netCDF4 as nc4
import matplotlib.pyplot as plt
import cmocean.cm as cm
from matplotlib.patches import Rectangle
import argparse

def read_llc_faces_geometry_to_faces(L0_input_dir,llc):
    XC_faces = {}
    YC_faces = {}
    for i in range(1,6):
        if i<3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(L0_input_dir, 'tile00' + str(i) + '.mitgrid'), llc, 3*llc)
        if i==3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(L0_input_dir, 'tile00' + str(i) + '.mitgrid'), llc, llc)
        if i>3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(L0_input_dir, 'tile00' + str(i) + '.mitgrid'), 3*llc, llc)
        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face
    return(XC_faces, YC_faces)

def read_mask_reference_from_nc_dict(nc_dict_file,mask_name):

    ds = nc4.Dataset(nc_dict_file)
    grp = ds.groups[mask_name]
    faces = grp.variables['source_faces'][:]
    rows = grp.variables['source_rows'][:]
    cols = grp.variables['source_cols'][:]
    ds.close()

    return(faces,rows,cols)

def read_grid_geometry_from_nc(config_dir,model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    ds.close()
    return(XC, YC, Depth)

def create_plot(config_dir, L1_model_name, L1_XC, L1_YC, L1_Depth, mask_names, mask_arrays):

    mask_colors = ['red', 'orange', 'green', 'purple']

    fig = plt.figure(figsize=(8, 6))
    plt.style.use('dark_background')

    # plt.subplot(1, 2, 1)
    plt.contour(L1_XC, L1_YC, L1_Depth, levels=[0, 500, 1000, 1500, 2000], colors='k', linewidths=0.25)
    C = plt.pcolormesh(L1_XC, L1_YC, L1_Depth, cmap=cm.deep, shading='nearest')  # ,vmin=vmin,vmax=vmax)
    plt.colorbar(C)
    plt.title(L1_model_name)
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    for m in range(len(mask_names)):
        plt.plot(mask_arrays[m][:,0], mask_arrays[m][:,1], '.', markersize=4, color=mask_colors[m])

    output_file = os.path.join(config_dir, 'L0', 'plots',
                               L1_model_name + '_dv_mask_locations.png')

    plt.savefig(output_file, bbox_inches='tight')
    plt.close(fig)

########################################################################################################################

def plot_dv_masks(config_dir,ecco_path,mask_names):
    llc = 270
    L1_model_name = 'L1_CE_Greenland'

    L0_input_dir = os.path.join(ecco_path,'LLC'+str(llc)+'_Files','mitgrid_tiles')
    XC_faces, YC_faces = read_llc_faces_geometry_to_faces(L0_input_dir,llc)

    nc_dict_file = os.path.join(config_dir,'L0','input','L0_dv_mask_reference_dict.nc')
    mask_arrays = []
    for mask in mask_names:
        mask_name = mask.split('_')[1]
        faces, rows, cols = read_mask_reference_from_nc_dict(nc_dict_file, mask_name)
        mask_points = np.zeros((len(faces),2))
        zero_ref_points = 0
        for f in range(len(faces)):
            if faces[f]!=0:
                mask_points[f, 0] = XC_faces[faces[f]][rows[f],cols[f]]
                mask_points[f, 1] = YC_faces[faces[f]][rows[f], cols[f]]
            else:
                print(mask,faces[f],rows[f],cols[f])
                zero_ref_points+=1
        if zero_ref_points>0:
            print(mask+' has some zero\'d reference points...')
        mask_points = mask_points[mask_points[:,0]!=0,:]
        mask_arrays.append(mask_points)

    L1_XC, L1_YC, L1_Depth = read_grid_geometry_from_nc(config_dir, L1_model_name)

    create_plot(config_dir, L1_model_name, L1_XC, L1_YC, L1_Depth, mask_names, mask_arrays)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The name of the configuration.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=True)

    parser.add_argument("-m", "--masks", action="store", help="List of masks to plot. "
                                                              "Default value is east south west.", default='', dest="masks", type=str, nargs='+',
                        required=False)

    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_path = args.ecco_path
    masks = args.masks

    if masks=='':
        masks = ['L1_west_BC_mask','L1_south_BC_mask','L1_east_BC_mask']

    plot_dv_masks(config_dir,ecco_path,masks)
