
import os
import numpy as np
import netCDF4 as nc4
from pyproj import Transformer
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import ast


def create_surface_hFacC_grid(bathy, delR, hFacMin, hFacMinDr):
    # This is from MITgcm inside ini_masks_etc.F
    # Note: R_low is just the bathy
    # Note: drF is the grid spacing (provided)
    # Note: recip_drF is the reciprocal of drF
    # Note: Ro_surf is the z coord of the surface (essentially always 0 for the ocean?)
    # C--   Calculate lopping factor hFacC : over-estimate the part inside of the domain
    # C     taking into account the lower_R Boundary (Bathymetry / Top of Atmos)
    #         DO k=1, Nr
    #          hFacMnSz = MAX( hFacMin, MIN(hFacMinDr*recip_drF(k),oneRL) )
    #          DO j=1-OLy,sNy+OLy
    #           DO i=1-OLx,sNx+OLx
    # C      o Non-dimensional distance between grid bound. and domain lower_R bound.
    #            hFac_loc = (rF(k)-R_low(i,j,bi,bj))*recip_drF(k)
    # C      o Select between, closed, open or partial (0,1,0-1)
    #            hFac_loc = MIN( MAX( hFac_loc, zeroRL ) , oneRL )
    # C      o Impose minimum fraction and/or size (dimensional)
    #            IF ( hFac_loc.LT.hFacMnSz*halfRL .OR.
    #      &          R_low(i,j,bi,bj).GE.Ro_surf(i,j,bi,bj) ) THEN
    #              hFacC(i,j,k,bi,bj) = zeroRS
    #            ELSE
    #              hFacC(i,j,k,bi,bj) = MAX( hFac_loc, hFacMnSz )
    #            ENDIF
    #           ENDDO
    #          ENDDO
    #         ENDDO

    # Define grids with same names as those in MITgcm
    RU = -1 * np.cumsum(delR)
    RL = np.concatenate([[0], RU[:-1]])
    R_low = bathy
    drF = RL - RU
    recip_drF = 1 / drF

    # Pythonize the above loops
    hFacC = np.zeros((np.shape(bathy)[0], np.shape(bathy)[1]))
    # for k in range(len(RL)):
    # if k%5==0:
    #     print('     - Calculating hFacC for depth cells '+str(k)+' to '+str(k+5))
    k=0
    hFacMnSz = np.max([hFacMin, np.min([hFacMinDr * recip_drF[k], 1])])
    for i in range(np.shape(bathy)[0]):
        for j in range(np.shape(bathy)[1]):
            #      o Non-dimensional distance between grid bound. and domain lower_R bound.
            hFac_loc = (RL[k] - R_low[i, j]) * recip_drF[k]
            #      o Select between, closed, open or partial (0,1,0-1)
            hFac_loc = np.min([np.max([hFac_loc, 0]), 1])
            #      o Impose minimum fraction and/or size (dimensional)
            if hFac_loc <= hFacMnSz * 0.5 or R_low[i, j] >= 0:
                hFacC[i, j] = 0
            else:
                hFacC[i, j] = np.max([hFac_loc, hFacMnSz])

    return(hFacC)

def generate_connected_mask(start_row, start_col, wet_grid):

    if wet_grid[start_row,start_col]==0:
        raise ValueError(' The start row/col location is  dry')

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

def create_bathymetry_file(config_dir, level_name, model_name,
                           central_wet_row, central_wet_col, hFacMinDr, hFacMin, delR, print_level=1):

    if print_level>=1:
        print('    - Generating the bathymetry for the ' + model_name + ' model from GEBCO and BedMachine')

    f = open(os.path.join(config_dir, 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    L_size = size_dict[model_name]
    n_rows = L_size[0]
    n_cols = L_size[1]

    file_path = os.path.join(config_dir, 'mitgrids', model_name + '.mitgrid')
    entire_grid = np.fromfile(file_path, dtype='>f8')
    entire_grid = np.reshape(entire_grid, (16, n_rows + 1, n_cols + 1))

    Lon_C = entire_grid[0, :, :]
    Lat_C = entire_grid[1, :, :]

    Lon_C = Lon_C[:-1, :-1]
    Lat_C = Lat_C[:-1, :-1]

    # x2, y2 = transformer.transform(polygon_array[:, y_column], polygon_array[:, x_column])
    transformer = Transformer.from_crs('EPSG:' + str(4326), 'EPSG:' + str(3413))
    XC, YC = transformer.transform(Lat_C.ravel(), Lon_C.ravel())

    XC = XC.reshape(np.shape(Lon_C))
    YC = YC.reshape(np.shape(Lat_C))

    ###############################################################
    # First interpolate with BedMachine

    if print_level >= 2:
        print('        - Interpolating from the BedMachine dataset')

    bathy = -10000*np.ones_like(XC)

    bedMachine_file = '/Users/michwood/Documents/Research/Data Repository/Greenland/Bathymetry/' \
                      'Bed Machine/Versions/BedMachineGreenland-2021-04-20.nc'

    ds = nc4.Dataset(bedMachine_file)
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    bed = ds.variables['bed'][:,:]
    surface = ds.variables['surface'][:,:]
    ds.close()

    min_x_index = np.argmin(np.abs(x - np.min(XC)))
    max_x_index = np.argmin(np.abs(x - np.max(XC)))
    max_y_index = np.argmin(np.abs(y - np.min(YC)))
    min_y_index = np.argmin(np.abs(y - np.max(YC)))

    x = x[min_x_index:max_x_index]
    y = y[min_y_index:max_y_index]
    bed = bed[min_y_index:max_y_index, min_x_index:max_x_index]
    surface = surface[min_y_index:max_y_index, min_x_index:max_x_index]

    # this is a manual nearest neighbor method
    bedmachine_resolution = 150
    for i in range(np.shape(XC)[0]):
        for j in range(np.shape(YC)[1]):
            x_index = np.argmin(np.abs(x - XC[i, j]))
            y_index = np.argmin(np.abs(y - YC[i, j]))
            dist = ((x[x_index]-XC[i,j])**2 + (y[y_index]-YC[i,j])**2)**0.5
            if dist<=bedmachine_resolution*np.sqrt(2)/2:
                if surface[y_index,x_index]>0:
                    bathy[i,j] = 0
                else:
                    bathy[i, j] = bed[y_index, x_index]

    # # this uses scipy
    # X,Y = np.meshgrid(x,y)
    # bed_interp = griddata(np.column_stack([X.ravel(),Y.ravel()]),bed.ravel(),(XC,YC))
    # bathy[~np.isnan(bed_interp)] = bed[~np.isnan(bed_interp)]

    # C = plt.imshow(bathy,origin='lower')
    # plt.colorbar(C)
    # plt.show()

    ###############################################################
    # Then interpolate with GEBCO if necessary

    if np.any(bathy<-9999):
        if print_level >= 2:
            print('        - Interpolating from the GEBCO dataset (not all points covered by BedMachine)')
        gebco_file = '/Users/michwood/Documents/Research/Data Repository/Global/Bathymetry/GEBCO/gebco_2021/GEBCO_2021.nc'
        ds = nc4.Dataset(gebco_file)
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        elev = ds.variables['elevation'][:,:]
        ds.close()

        min_lon_index = np.argmin(np.abs(lon - np.min(Lon_C)))
        max_lon_index = np.argmin(np.abs(lon - np.max(Lon_C)))
        min_lat_index = np.argmin(np.abs(lat - np.min(Lat_C)))
        max_lat_index = np.argmin(np.abs(lat - np.max(Lat_C)))

        lon = lon[min_lon_index:max_lon_index]
        lat = lat[min_lat_index:max_lat_index]
        elev = elev[min_lat_index:max_lat_index, min_lon_index:max_lon_index]

        # plt.imshow(elev,origin='lower')
        # plt.show()

        # this is a manual nearest neighbor method
        for i in range(np.shape(XC)[0]):
            for j in range(np.shape(YC)[1]):
                lon_index = np.argmin(np.abs(lon - Lon_C[i, j]))
                lat_index = np.argmin(np.abs(lat - Lat_C[i, j]))
                if elev[lat_index,lon_index]<0:
                    bathy[i,j] = elev[lat_index,lon_index]
                else:
                    bathy[i,j] = 0

    C = plt.imshow(bathy, origin='lower')
    plt.colorbar(C)
    plt.title(np.shape(bathy))
    plt.show()

    ###############################################################
    # Next, make some adjustments to the bathymetry

    min_cell_height = 10
    if print_level >= 2:
        print('        - Lowering shallow points down to a depth of '+str(min_cell_height)+' m')
    bathy[np.logical_and(bathy>-min_cell_height,bathy<0)] = -min_cell_height

    # C = plt.imshow(bathy, origin='lower')
    # plt.colorbar(C)
    # plt.title(np.shape(bathy))
    # plt.show()

    if print_level >= 2:
        print('        - Generating a wet grid for the domain')
    hFacC_grid = create_surface_hFacC_grid(bathy, delR, hFacMin, hFacMinDr)

    wet_grid = np.copy(hFacC_grid)
    wet_grid[wet_grid>1] = 1

    # C = plt.imshow(wet_grid, origin='lower')
    # plt.colorbar(C)
    # plt.title(np.shape(bathy))
    # plt.show()

    if print_level >= 2:
        print('        - Generating a mask for the domain to eliminate unconnected regions')
    mask_grid = generate_connected_mask(central_wet_row, central_wet_col, wet_grid)

    # C = plt.imshow(mask_grid, origin='lower')
    # plt.colorbar(C)
    # plt.title(np.shape(bathy))
    # plt.show()

    bathy[mask_grid==0] = 0

    # C = plt.imshow(bathy, origin='lower')
    # plt.colorbar(C)
    # plt.title(np.shape(bathy))
    # plt.show()

    output_dir = os.path.join(config_dir,level_name,model_name)
    if 'input' not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,'input'))

    output_dir = os.path.join(output_dir,'input')

    output_file = os.path.join(output_dir, model_name+'_bathymetry.bin')
    bathy.ravel('C').astype('>f4').tofile(output_file)

