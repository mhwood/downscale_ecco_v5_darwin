import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from math import radians, cos, sin, asin, sqrt
import simplegrid as sg
from scipy.interpolate import griddata
from ecco_v4_py.llc_array_conversion import llc_faces_to_compact, llc_compact_to_faces
import argparse
import ast


def read_global_XC_YC_bathy_HFac(ecco_dir,llc):

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

    hFac_compact = np.fromfile(os.path.join(input_init_dir, 'hFacC.data'), '>f4')
    hFac_compact = np.reshape(hFac_compact, (50, 13 * llc, llc))
    hFac_faces = llc_compact_to_faces(hFac_compact, less_output=True)

    return(XC_faces,YC_faces, bathy_faces, hFac_faces)

def read_mask_compact_to_faces(config_dir, file_name, llc):
    mask_compact = np.fromfile(os.path.join(config_dir,'L0','input','dv',file_name), '>f4')
    mask_compact = np.reshape(mask_compact, (13 * llc, llc))
    mask_faces = llc_compact_to_faces(mask_compact, less_output=True)
    return(mask_faces)

def count_mask_points_per_proc(mask_faces,bathy_faces,hFac_faces,sNx):

    proc_counter = 0
    max_points_per_proc = 0
    max_depth_cell_for_mask = 0
    for face in range(1,6):
        mask_face = mask_faces[face]
        bathy_face = bathy_faces[face]
        hFac_face = hFac_faces[face]
        ll_row = 0
        ll_col = 0
        procs_in_face = int(np.size(mask_face)/(sNx*sNx))
        for proc in range(procs_in_face):
            mask_proc = mask_face[ll_row:ll_row+sNx,ll_col:ll_col+sNx]
            bathy_proc = bathy_face[ll_row:ll_row + sNx, ll_col:ll_col + sNx]
            hFac_proc = hFac_face[:,ll_row:ll_row + sNx, ll_col:ll_col + sNx]
            # if np.any(bathy_proc!=0):
            proc_counter+=1
            if np.any(mask_proc>0):
                mask_proc_points = np.sum(mask_proc!=0)
                if max_points_per_proc<mask_proc_points:
                    max_points_per_proc = mask_proc_points
                # print('            - Proc '+str(proc_counter)+' has '+str(mask_proc_points)+' points')
                nz_rows,nz_cols = np.where(mask_proc!=0)
                for ri in range(len(nz_rows)):
                    hFac = hFac_proc[:,nz_rows[ri],nz_cols[ri]]
                    nonzero_levels = np.sum(hFac!=0)
                    if nonzero_levels>max_depth_cell_for_mask:
                        max_depth_cell_for_mask = nonzero_levels
            if ll_col+sNx>=np.shape(mask_face)[1]:
                ll_row+=sNx
                ll_col=0
            else:
                ll_col+=sNx

    # print('Total procs:',proc_counter)
    print('          - Max points per proc: '+str(max_points_per_proc))
    print('          - Deepest depth layer: ' + str(max_depth_cell_for_mask))


    return(max_points_per_proc, max_depth_cell_for_mask)

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

def check_dv_mask_proc_points(config_dir, ecco_path,L1_model_names,print_status):

    if print_status:
        print('Checking the number of points per processor in each diagnostics_vec mask')

    llc = 270
    sNx = 30

    if print_status:
        print('    - Reading in the L0 domain files')

    # read the mitgrids to faces
    if print_status:
        print('    - Reading the LLC '+str(llc)+' geometry')
    ecco_dir = os.path.join(ecco_path,'LLC'+str(llc)+'_Files')
    XC_faces, YC_faces, bathy_faces, hFac_faces = read_global_XC_YC_bathy_HFac(ecco_dir, llc)

    for L1_model_name in L1_model_names:
        print('    - Checking masks for the '+L1_model_name+' model')

        for boundary in ['west','north','south','east']:
            file_name = L1_model_name+'_'+boundary+'.bin'
            if file_name in os.listdir(os.path.join(config_dir,'L0','input','dv')):
                print('        - Reading file '+file_name)
                mask_faces = read_mask_compact_to_faces(config_dir, file_name, llc)

                max_points_per_proc, max_depth_cell_for_mask = \
                    count_mask_points_per_proc(mask_faces, bathy_faces, hFac_faces, sNx)











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

    check_dv_mask_proc_points(config_dir,ecco_path,L1_model_names,print_status)
