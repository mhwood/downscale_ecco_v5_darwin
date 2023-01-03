
import os
import shutil
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import argparse
import ast


def read_grid_information(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    rA = np.array(ds.variables['rA'][:, :])
    dxG = np.array(ds.variables['dxG'][:, :])
    dyG = np.array(ds.variables['dyG'][:, :])
    hFacS = np.array(ds.variables['HFacS'][:, :, :])
    hFacW = np.array(ds.variables['HFacW'][:, :, :])
    delR = np.array(ds.variables['drF'][:])
    ds.close()
    return (rA, dxG, dyG, hFacS, hFacW, delR)


def read_velocities_normal_to_boundary(config_dir, model_name,year,boundary,Nr,n_cols,n_rows,balanced):

    if balanced:
        balanced_text = ''
    else:
        balanced_text = '_unbalanced'

    obcs_folder = os.path.join(config_dir,'L1', model_name,'input', 'obcs')

    if boundary=='east':

        uvel_east_file = os.path.join(obcs_folder,'L1_BC_east_UVEL_'+str(year)+balanced_text)
        uvel_east = np.fromfile(uvel_east_file,'>f4')
        n_timesteps = int(np.size(uvel_east)/(Nr*n_rows))
        uvel_east = np.reshape(uvel_east,(n_timesteps,Nr,n_rows)).astype(float)
        vel_out = uvel_east

    if boundary=='west':

        uvel_west_file = os.path.join(obcs_folder,'L1_BC_west_UVEL_'+str(year)+balanced_text)
        uvel_west = np.fromfile(uvel_west_file,'>f4')
        n_timesteps = int(np.size(uvel_west)/(Nr*n_rows))
        uvel_west = np.reshape(uvel_west,(n_timesteps,Nr,n_rows)).astype(float)
        vel_out = uvel_west

    # plt.subplot(1,3,1)
    # C = plt.imshow(uvel_east[17,:,:])
    # plt.colorbar(C)
    # plt.subplot(1, 3, 2)
    # C = plt.imshow(uvel_east[23, :, :])
    # plt.colorbar(C)
    # plt.subplot(1, 3, 3)
    # C = plt.imshow(uvel_east[23, :, :]!=0)
    # plt.colorbar(C)
    # plt.show()

    if boundary == 'north':

        vvel_north_file = os.path.join(obcs_folder, 'L1_BC_north_VVEL_'+str(year)+balanced_text)
        vvel_north = np.fromfile(vvel_north_file, '>f4')
        n_timesteps = int(np.size(vvel_north) / (Nr * n_cols))
        vvel_north = np.reshape(vvel_north, (n_timesteps, Nr, n_cols)).astype(float)
        vel_out = vvel_north

    if boundary == 'south':

        vvel_south_file = os.path.join(obcs_folder, 'L1_BC_south_VVEL_'+str(year)+balanced_text)
        vvel_south = np.fromfile(vvel_south_file, '>f4')
        n_timesteps = int(np.size(vvel_south) / (Nr * n_cols))
        vvel_south = np.reshape(vvel_south, (n_timesteps, Nr, n_cols)).astype(float)
        vel_out = vvel_south

    return(vel_out)


def calculate_flux(boundary,timestep,field,width,delR,hFac,print_snapshot_stats):
    flux = np.zeros_like(field).astype(float)
    total_area = 0.0
    for i in range(np.shape(field)[1]):
        for j in range(np.shape(field)[0]):
            flux[j,i] = width[i]*delR[j]*hFac[j,i]*field[j,i]
            total_area += width[i]*delR[j]*hFac[j,i]
    total_flux = np.sum(flux)
    return(total_flux,total_area)


def calculate_flux_timeseries(vel_grid, boundary,
                              dxG, dyG, hFacS, hFacW, delR, print_level, print_snapshot_stats = False):

    if boundary=='north':
        hFac = hFacS[:,-2,:]
        width = dxG[-2,:]
    if boundary=='south':
        hFac = hFacS[:,1,:]
        width = dxG[1, :]
    if boundary=='east':
        hFac = hFacW[:,:,-2]
        width = dyG[:, -2]
    if boundary=='west':
        hFac = hFacW[:,:,1]
        width = dyG[:, 1]

    # note: the first and last grid cells are not calculated in the flux
    # one way to solve this issue is just to make hfac=0 in these areas
    hFac[:, 0] = 0
    hFac[:, -1] = 0

    # make a mask to apply the corrections later
    mask = np.copy(hFac)
    mask[mask>0]=1
    mask[mask<=0]=0

    n_timesteps = np.shape(vel_grid)[0]
    flux_timeseries = np.zeros((n_timesteps,)).astype(float)

    for timestep in range(n_timesteps):
        if print_level>=4:
            if timestep%1000==0:
                print('                - Working on timesteps '+str(timestep)+' to '+str(np.min([timestep+1000,n_timesteps]))+
                      ' out of '+str(n_timesteps)+' for the '+boundary+' boundary')
        total_flux, total_area = calculate_flux(boundary, timestep, vel_grid[timestep,:,:], width, delR, hFac, print_snapshot_stats)
        flux_timeseries[timestep] = total_flux

    if print_snapshot_stats:
        print('    Area ',total_area)
        print('    flow ',flux_timeseries[0])

    if boundary=='north':
        flux_timeseries *= -1
    if boundary=='east':
        flux_timeseries *= -1

    return(flux_timeseries, total_area, mask)


def plot_boundary_fluxes(output_file, boundaries, east_flux_timeseries, west_flux_timeseries, north_flux_timeseries,south_flux_timeseries,rA):
    flux_started = False
    if 'east' in boundaries:
        if not flux_started:
            flux_started = True
            total_flux_timeseries = east_flux_timeseries
        else:
            total_flux_timeseries = total_flux_timeseries + east_flux_timeseries
    if 'west' in boundaries:
        if not flux_started:
            flux_started = True
            total_flux_timeseries = west_flux_timeseries
        else:
            total_flux_timeseries = total_flux_timeseries + west_flux_timeseries
    if 'north' in boundaries:
        if not flux_started:
            flux_started = True
            total_flux_timeseries = north_flux_timeseries
        else:
            total_flux_timeseries = total_flux_timeseries + north_flux_timeseries
    if 'south' in boundaries:
        if not flux_started:
            flux_started = True
            total_flux_timeseries = south_flux_timeseries
        else:
            total_flux_timeseries = total_flux_timeseries + south_flux_timeseries

    integrated_flux = np.cumsum(3600*total_flux_timeseries)

    fig = plt.figure(figsize=(12, 14))
    plt.style.use("dark_background")

    for b in range(len(boundaries)):
        boundary = boundaries[b]

        plt.subplot(len(boundaries)+2, 1, b+1)
        if boundary == 'east':
            plt.plot(east_flux_timeseries/1e6)
            plt.title('East Boundary Flux Into Domain')
        if boundary == 'west':
            plt.plot(west_flux_timeseries/1e6)
            plt.title('West Boundary Flux Into Domain')
        if boundary == 'south':
            plt.plot(south_flux_timeseries/1e6)
            plt.title('South Boundary Flux Into Domain')
        if boundary == 'north':
            plt.plot(north_flux_timeseries/1e6)
            plt.title('North Boundary Flux Into Domain')
        plt.grid(linestyle='--', alpha=0.2)
        plt.gca().set_xticklabels([])
        plt.ylabel('Sv')

    plt.subplot(len(boundaries)+2, 1, len(boundaries)+1)
    plt.plot(total_flux_timeseries/1e6)
    plt.grid(linestyle='--', alpha=0.2)
    plt.title('Total Flux Into Domain')
    plt.gca().set_xticklabels([])
    plt.ylabel('Sv')

    plt.subplot(len(boundaries)+2, 1, len(boundaries)+2)
    plt.plot(integrated_flux/np.sum(rA))
    plt.title('Integrated Flux Into Domain')
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('m (averaged over domain)')

    plt.savefig(output_file,bbox_inches='tight')
    plt.close(fig)


def balance_flux_on_boundaries(config_dir, boundaries, model_name, year, running_average_radius, Nr,n_cols,n_rows,
                               east_flux_timeseries, east_area, east_mask,
                               west_flux_timeseries, west_area, west_mask,
                               north_flux_timeseries, north_area, north_mask,
                               south_flux_timeseries, south_area, south_mask,
                               print_level):

    #####################################################
    # calculate the flux adjustment on each boundary
    # following the obcs package, a constant value is removed from all wet cells on the boundary
    # one can modify the balance fractions below, but the default it to remove it from all boundaries equally
    # note, for these, that the W and S boundaries have a -1 applied below

    if 'north' in boundaries:
        balanceFacN = 1
    if 'south' in boundaries:
        balanceFacS = 1
    if 'east' in boundaries:
        balanceFacE = 1
    if 'west' in boundaries:
        balanceFacW = 1

    flux_started = False
    if 'east' in boundaries:
        if not flux_started:
            flux_started = True
            total_flux_timeseries = east_flux_timeseries
        else:
            total_flux_timeseries = total_flux_timeseries + east_flux_timeseries
    if 'west' in boundaries:
        if not flux_started:
            flux_started = True
            total_flux_timeseries = west_flux_timeseries
        else:
            total_flux_timeseries = total_flux_timeseries + west_flux_timeseries
    if 'north' in boundaries:
        if not flux_started:
            flux_started = True
            total_flux_timeseries = north_flux_timeseries
        else:
            total_flux_timeseries = total_flux_timeseries + north_flux_timeseries
    if 'south' in boundaries:
        if not flux_started:
            flux_started = True
            total_flux_timeseries = south_flux_timeseries
        else:
            total_flux_timeseries = total_flux_timeseries + south_flux_timeseries

    total_area = east_area + north_area + south_area + west_area
    total_correction = total_flux_timeseries / total_area

    print('|- corrections -----------------------------|')
    if 'east' in boundaries:
        print('    flowE ',east_flux_timeseries[0])
    if 'west' in boundaries:
        print('    flowE ',west_flux_timeseries[0])
    if 'north' in boundaries:
        print('    flowN ', north_flux_timeseries[0])
    if 'south' in boundaries:
        print('    flowS ', south_flux_timeseries[0])
    print('    inFlow ',total_flux_timeseries[0])
    print('    areaOB ',total_area)
    print('    inFlow/areaOB ',total_correction[0])

    if running_average_radius>0:
        total_correction_copy = np.copy(total_correction)
        smooth_total_correction = np.zeros((len(total_correction),))
        for i in range(len(total_correction)):
            if i>=running_average_radius:
                smooth_total_correction[i] = np.mean(total_correction[i-running_average_radius:i])
        smooth_total_correction[:running_average_radius] = smooth_total_correction[running_average_radius]

        total_correction = smooth_total_correction

        # plt.plot(total_correction_copy,'b-')
        # plt.plot(smooth_total_correction,'g-')
        # plt.show()

    #####################################################
    # apply the boundary flux adjustments to each boundary

    n_timesteps = np.size(total_flux_timeseries)

    if 'east' in boundaries:
        print('    - Correcting the east flux')
        uvel_east = read_velocities_normal_to_boundary(config_dir,model_name,  year,'east', Nr, n_cols, n_rows, balanced=False)
        for i in range(n_timesteps):
            uvel_east[i,:,:] += total_correction[i] * east_mask*balanceFacE
        output_file = os.path.join(config_dir,'L1',model_name, 'input', 'obcs', 'L1_BC_east_UVEL_'+str(year))
        uvel_east.ravel('C').astype('>f4').tofile(output_file)
        del uvel_east

    if 'west' in boundaries:
        print('    - Correcting the west flux')
        uvel_west = read_velocities_normal_to_boundary(config_dir,model_name,  year,'west', Nr, n_cols, n_rows, balanced=False)
        for i in range(n_timesteps):
            uvel_west[i,:,:] += total_correction[i] * west_mask*balanceFacW
        output_file = os.path.join(config_dir,'L1',model_name, 'input', 'obcs', 'L1_BC_west_UVEL_'+str(year))
        uvel_west.ravel('C').astype('>f4').tofile(output_file)
        del uvel_west

    if 'north' in boundaries:
        print('    - Correcting the north flux')
        vvel_north = read_velocities_normal_to_boundary(config_dir,model_name,  year,'north', Nr, n_cols, n_rows, balanced=False)
        for i in range(n_timesteps):
            vvel_north[i, :, :] += total_correction[i] * north_mask * balanceFacN
        output_file = os.path.join(config_dir,'L1',model_name, 'input', 'obcs', 'L1_BC_north_VVEL_'+str(year))
        vvel_north.ravel('C').astype('>f4').tofile(output_file)
        del vvel_north

    if 'south' in boundaries:
        print('    - Correcting the south flux')
        vvel_south = read_velocities_normal_to_boundary(config_dir,model_name,  year,'south', Nr, n_cols, n_rows, balanced=False)
        for i in range(n_timesteps):
            vvel_south[i, :, :] -= total_correction[i] * south_mask * balanceFacS
        output_file = os.path.join(config_dir,'L1',model_name, 'input', 'obcs', 'L1_BC_south_VVEL_'+str(year))
        vvel_south.ravel('C').astype('>f4').tofile(output_file)
        del vvel_south

def balance_bc_fields(config_dir, model_name, boundaries, start_year, final_year, print_level):

    if print_level>=1:
        print('    - Balance the fluxes into the model domain to correct for interpolation biases')

    if print_level >= 2:
        print('         - Reading in the grid information')
    rA, dxG, dyG, hFacS, hFacW, delR = read_grid_information(config_dir, model_name)
    n_rows = np.shape(rA)[0]
    n_cols = np.shape(rA)[1]
    Nr = len(delR)

    running_average_radius = 0

    for year in range(start_year,final_year+1):

        if print_level >= 1:
            print('    - Working on year '+str(year))

        ########################################################################################################

        if 'east' in boundaries:

            if print_level >= 2:
                print('        - Working on the east boundary')

            if print_level >= 3:
                print('            - Reading in the east velocity grid')
            uvel_east = read_velocities_normal_to_boundary(config_dir, model_name, year,'east',Nr,n_cols,n_rows,balanced=False)

            if print_level >= 3:
                print('            - Calculating the fluxes on the east boundary')
            east_flux_timeseries, east_area, east_mask  = \
                calculate_flux_timeseries(uvel_east, 'east', dxG, dyG, hFacS, hFacW, delR, print_level, print_snapshot_stats = False)

            del uvel_east

        else:
            east_flux_timeseries = []
            east_area = 0
            east_mask = 0

        ########################################################################################################

        if 'west' in boundaries:

            if print_level >= 2:
                print('        - Working on the west boundary')

            if print_level >= 3:
                print('            - Reading in the west velocity grid')
            uvel_west = read_velocities_normal_to_boundary(config_dir, model_name, year, 'west', Nr, n_cols, n_rows,
                                                           balanced=False)

            if print_level >= 3:
                print('            - Calculating the fluxes on the west boundary')
            west_flux_timeseries, west_area, west_mask = \
                calculate_flux_timeseries(uvel_west, 'west', dxG, dyG, hFacS, hFacW, delR, print_level,
                                          print_snapshot_stats=False)

            del uvel_west

        else:
            west_flux_timeseries = []
            west_area = 0
            west_mask = 0

        ########################################################################################################

        if 'north' in boundaries:
            if print_level >= 2:
                print('         - Working on the north boundary')

            if print_level >= 3:
                print('            - Reading in the north velocity grid')
            vvel_north = read_velocities_normal_to_boundary(config_dir, model_name,year,'north', Nr, n_cols, n_rows,balanced=False)

            if print_level >= 3:
                print('            - Calculating the fluxes on the north boundary')
            north_flux_timeseries, north_area, north_mask = \
                calculate_flux_timeseries(vvel_north, 'north', dxG, dyG, hFacS, hFacW, delR, print_level, print_snapshot_stats=False)

            del vvel_north

        else:
            north_flux_timeseries = []
            north_area = 0
            north_mask = 0

        ########################################################################################################

        if 'south' in boundaries:
            if print_level >= 2:
                print('         - Working on the south boundary')

            if print_level >= 3:
                print('            - Reading in the south velocity grid')
            vvel_south = read_velocities_normal_to_boundary(config_dir, model_name, year,'south', Nr, n_cols, n_rows,balanced=False)

            if print_level >= 3:
                print('            - Calculating the fluxes on the south boundary')
            south_flux_timeseries, south_area, south_mask = \
                calculate_flux_timeseries(vvel_south, 'south', dxG, dyG, hFacS, hFacW, delR, print_level, print_snapshot_stats=False)

            del vvel_south

        else:
            north_flux_timeseries = []
            north_area = 0
            north_mask = 0

        ########################################################################################################

        print(' - Step 5: Plotting the uncorrected fluxes')
        output_file = os.path.join(config_dir,'L1',model_name,'plots','init_files','BCs',model_name+'_unbalanced_fluxes_'+str(year)+'.png')
        plot_boundary_fluxes(output_file, boundaries, east_flux_timeseries, west_flux_timeseries, north_flux_timeseries, south_flux_timeseries,rA)

        print(' - Step 6: Adjusting the fluxes on each boundary')
        balance_flux_on_boundaries(config_dir, boundaries, model_name, year, running_average_radius, Nr,n_cols,n_rows,
                                   east_flux_timeseries, east_area, east_mask,
                                   west_flux_timeseries, west_area, west_mask,
                                   north_flux_timeseries, north_area, north_mask,
                                   south_flux_timeseries, south_area, south_mask,
                                   print_level)

        # ########################################################################################################
        #
        # if print_level >= 2:
        #     print('        - Working on the balanced east boundary')
        #
        # if print_level >= 3:
        #     print('            - Reading in the east velocity grid')
        # uvel_east = read_velocities_normal_to_boundary(config_dir, model_name, year, 'east', Nr, n_cols, n_rows,
        #                                                balanced=True)
        #
        # if print_level >= 3:
        #     print('            - Calculating the fluxes on the east boundary')
        # east_flux_timeseries, east_area, east_mask = \
        #     calculate_flux_timeseries(uvel_east, 'east', dxG, dyG, hFacS, hFacW, delR, print_level,
        #                               print_snapshot_stats=False)
        #
        # del uvel_east
        #
        # ########################################################################################################
        #
        # if print_level >= 2:
        #     print('         - Working on the balanced north boundary')
        #
        # if print_level >= 3:
        #     print('            - Reading in the north velocity grid')
        # vvel_north = read_velocities_normal_to_boundary(config_dir, model_name, year, 'north', Nr, n_cols, n_rows,
        #                                                 balanced=True)
        #
        # if print_level >= 3:
        #     print('            - Calculating the fluxes on the north boundary')
        # north_flux_timeseries, north_area, north_mask = \
        #     calculate_flux_timeseries(vvel_north, 'north', dxG, dyG, hFacS, hFacW, delR, print_level,
        #                               print_snapshot_stats=False)
        #
        # del vvel_north
        #
        # ########################################################################################################
        #
        # if print_level >= 2:
        #     print('         - Working on the balanced south boundary')
        #
        # if print_level >= 3:
        #     print('            - Reading in the south velocity grid')
        # vvel_south = read_velocities_normal_to_boundary(config_dir, model_name, year, 'south', Nr, n_cols, n_rows,
        #                                                 balanced=True)
        #
        # if print_level >= 3:
        #     print('            - Calculating the fluxes on the south boundary')
        # south_flux_timeseries, south_area, south_mask = \
        #     calculate_flux_timeseries(vvel_south, 'south', dxG, dyG, hFacS, hFacW, delR, print_level,
        #                               print_snapshot_stats=False)
        #
        # del vvel_south
        #
        # ########################################################################################################
        #
        # print(' - Step 6: Plotting the corrected fluxes')
        # output_file = os.path.join(config_dir,'L1',model_name,'plots','analysis','flux_balances',model_name+'_balanced_fluxes_'+str(year)+'.png')
        # plot_boundary_fluxes(output_file, east_flux_timeseries,north_flux_timeseries,south_flux_timeseries,rA)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-m", "--L1_model_name", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="L1_model_name",
                        type=str, required=True)

    parser.add_argument("-S", "--start_year", action="store",
                        help="The start year.", dest="start_year",
                        type=int, required=True)

    parser.add_argument("-F", "--final_year", action="store",
                        help="The final year.", dest="final_year",
                        type=int, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    L1_model_name = args.L1_model_name
    start_year = args.start_year
    final_year = args.final_year

    balance_bc_fields(config_dir, L1_model_name, start_year, final_year, print_level=4)
