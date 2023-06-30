
import os
import argparse
import sys
import numpy as np

def create_L1_W_Greenland_files(config_dir, ecco_dir):

    L1_model_name = 'L1_W_Greenland'

    parent_model_level = 'L0'
    parent_model_pickup_iteration = 73

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'init_file_creation'))
    sys.path.insert(1, os.path.join(config_dir, 'L1', L1_model_name, 'utils', 'init_file_creation'))

    print_level = 5

    steps = [3,4, 5, 6]
    steps = [2]

    # step 1: make the grid
    if 1 in steps:
        print('Step 1: Creating the mitgrid file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_mitgrid as cm
        cm.create_L1_W_Greenland_mitgrid_file(config_dir, L1_model_name, ecco_dir, print_level)

    # step 2: make the bathymetry
    if 2 in steps:
        print('Step 2: Creating the bathymetry file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_bathymetry as cb
        cb.create_L1_W_Greenland_bathymetry(config_dir, ecco_dir, print_level)

    #########################################################################################################
    # After the bathymetry is made, the *_for_grid model should be run to make the grid
    # This bit checks that the grid file has been made
    if np.any(np.array(steps)>2):
        grid_path = os.path.join(config_dir, 'nc_grids', L1_model_name+'_grid.nc')
        if not os.path.exists(grid_path):
            sys.exit("Need to make the grid netcdf file for this model for reference\n"
                     "   - Run the stitch_L1_W_Greenland_nc_grid_files_for_ref.py script with the -f 1 flag (and -r 150 -c 210)\n"
                     "   - Run the *for_grid model\n"
                     "   - Run the stitch_L1_W_Greenland_nc_grid_files_for_ref.py script")

    # step 3: make the initial conditions
    if 3 in steps:
        print('Step 3: Creating the pickup (initial conditions) file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_pickup as cp
        cp.create_pickup_file(config_dir, ecco_dir, print_level)

    # step 4: make the seaice initial conditions
    if 4 in steps:
        print('Step 4: Creating the seaice pickup (initial conditions) file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_seaice_pickup as cp
        cp.create_seaice_pickup_file(config_dir, ecco_dir, print_level)

    # step 5: make the seaice initial conditions
    if 5 in steps:
        print('Step 5: Creating the ptracer pickup (initial conditions) file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_ptracer_pickup as cp
        cp.create_ptracer_pickup_file(config_dir, ecco_dir, print_level)

    # step 6: make the darwin initial conditions
    if 6 in steps:
        print('Step 6: Creating the darwin pickup (initial conditions) file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_darwin_pickup as cp
        cp.create_darwin_pickup_file(config_dir, ecco_dir, print_level)

    #########################################################################################################

    # step 7: make the external forcing conditions
    if 7 in steps:
        print('Step 7: Creating the external forcing conditions for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_exf as ce
        ce.create_exf_files(config_dir, L1_model_name, parent_model_level, parent_model_name, print_level)

    # step 8: make the boundary conditions
    if 8 in steps:
        print('Step 8: Creating the boundary conditions for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_BCs as cbc
        cbc.create_BCs(config_dir, parent_model_name, L1_model_name, print_level)

    # step 9: make the dv masks
    if 9 in steps:
        print('Step 9: Creating the diagnostic_vec masks for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_dv_masks as cdv
        cdv.create_dv_masks(config_dir, L1_model_name, print_level)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_dir", action="store",
                        help="The directory where ECCO files are stored.", dest="ecco_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_dir = args.ecco_dir

    create_L1_W_Greenland_files(config_dir, ecco_dir)
   

