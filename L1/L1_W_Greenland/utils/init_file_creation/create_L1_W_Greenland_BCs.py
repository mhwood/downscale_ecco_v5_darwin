
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys
import ast


def create_BCs(config_dir, ecco_dir, print_level):

    L1_model_name = 'L1_W_Greenland'

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils','init_file_creation'))

    start_year = 1992
    start_month = 1

    final_year = 2018
    final_month = 12

    # # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    # import create_L1_BC_field_ref as ebcr
    # ebcr.create_L1_BC_ref_file(config_dir, L1_model_name, print_level)

    ##################################################################################################
    # In the first part of this script, we take an initial crack at the interpolation
    # generating monthly interpolated files for all variables

    proc_ids = np.arange(160).tolist() #(31+9)*4 = 160
    boundaries = ['north','south','west','east']

    # import create_L1_monthly_BCs as cef
    # for proc_id in proc_ids:
    #     cef.create_bc_fields_via_interpolation(config_dir, ecco_dir, L1_model_name, boundaries, proc_id,
    #                                            start_year, final_year, start_month, final_month, print_level)

    # #################################################################################################
    # # In the second part of this script, we balance the fluxes of velocity across the boundaries
    # # and ensure that the
    #
    # # step 3: combine all of the velocity BC fields into a single file
    # import combine_and_rotate_L1_monthly_bc_files as com
    # for proc_id in proc_ids:
    #     com.combine_and_rotate_L1_monthly_bcs(config_dir, L1_model_name, boundaries, proc_id,
    #                                         start_year, final_year, print_level, velocity_only = True)
    #
    # # step 4: generate the unbalanced timeseries for volume flux across the boundaries
    # import generate_L1_flux_timeseries as gL1t
    # for proc_id in [2]: # velocity only
    #     gL1t.create_L1_timeseries(config_dir, L1_model_name, boundaries, proc_id,
    #                               start_year, final_year, balanced=False, print_level=4)
    #
    # # step 5: calculate the velocity fluxes on the domain in the L0 model
    # import generate_L0_flux_timeseries_on_L1_boundary as gL0t
    # for proc_id in [2]: # velocity only
    #     gL0t.create_L0_timeseries(config_dir, ecco_dir, L1_model_name, boundaries, proc_id,
    #                                        start_year, final_year, print_level=4)
    #
    # # # step 6: remove anomalies and balance all of the velocity BC fields
    # # import balance_L1_vector_BCs as bbcs
    # # bbcs.balance_bc_fields(config_dir, L1_model_name, boundaries, start_year, final_year, print_level)
    #
    # # step 7: generate the timeseries for volume flux across the boundaries
    # import generate_L1_flux_timeseries as gL1t
    # for proc_id in [2]: # velocity only
    #     gL1t.create_L1_timeseries(config_dir, L1_model_name, boundaries, proc_id,
    #                               start_year, final_year, balanced=True, print_level=4)
    #
    # sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'plot_creation','init_files'))
    # import plot_L1_vs_L0_BC_flux_timeseries as pbc
    # pbc.plot_BC_timeseries_comparison(config_dir, L1_model_name, boundaries, velocity_only = True)

    # ##################################################################################################
    # # In the third part of this script, we balance the fluxes of all other tracers into the domain
    #
    # timeseries_proc_ids = list(np.arange(1,41))
    # timeseries_proc_ids.pop(timeseries_proc_ids.index(2)) #VEL
    # timeseries_proc_ids.pop(timeseries_proc_ids.index(3)) #UICE
    # timeseries_proc_ids.pop(timeseries_proc_ids.index(4)) #VICE
    # timeseries_proc_ids.pop(timeseries_proc_ids.index(5)) #HSNOW
    # timeseries_proc_ids.pop(timeseries_proc_ids.index(6)) #HEFF
    # timeseries_proc_ids.pop(timeseries_proc_ids.index(7)) #AREA
    #
    # # step 5: calculate the fluxes on the domain in the L0 model
    # import generate_L0_flux_timeseries_on_L1_boundary as gL0t
    # for proc_id in timeseries_proc_ids:
    #     gL0t.create_L0_timeseries(config_dir, ecco_dir, L1_model_name, boundaries, proc_id,
    #                                        start_year, final_year, print_level=4)
    #
    # # step 6: calculate the fluxes on the domain in the L1 model
    # import generate_L1_flux_timeseries as gL1t
    # for proc_id in timeseries_proc_ids:
    #     gL1t.create_L1_timeseries(config_dir, L1_model_name, boundaries, proc_id,
    #                               start_year, final_year, balanced = False, print_level=4)
    #
    # # step 7: balance boundary fluxes
    # import balance_L1_BCs_to_L0 as bal10
    # for proc_id in timeseries_proc_ids:
    #     bal10.balance_bc_fields(config_dir, L1_model_name, boundaries, proc_id,  start_year, final_year, print_level)
    #
    # # step 8: calculate the balanced fluxes on the domain in the L1 model
    # import generate_L1_flux_timeseries as gL1t
    # for proc_id in timeseries_proc_ids:
    #     gL1t.create_L1_timeseries(config_dir, L1_model_name, boundaries, proc_id,
    #                               start_year, final_year, balanced=True, print_level=4)

    # sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'plot_creation','init_files'))
    # import plot_L1_vs_L0_BC_flux_timeseries as pbc
    # pbc.plot_BC_timeseries_comparison(config_dir, L1_model_name, boundaries, velocity_only = False)

    ##################################################################################################
    # The last part of this script allows you to make annual files by combining without correcting

    import combine_and_rotate_L1_monthly_bc_files as cbc
    # for proc_id in [4,5,6,7,8]:
    var_names = ['UICE', 'VICE', 'HSNOW', 'HEFF', 'AREA']
    cbc.combine_and_rotate_L1_monthly_bcs(config_dir, L1_model_name, boundaries, var_names,
                                          start_year, final_year, print_level)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_dir", action="store",
                        help="The directory where ECCO files are stored.", dest="ecco_dir",
                        type=str, required=True)


    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_dir = args.ecco_dir

    create_BCs(config_dir, ecco_dir, print_level=4)
   

