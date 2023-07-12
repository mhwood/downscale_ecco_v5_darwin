import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import simplegrid as sg
from scipy.interpolate import griddata
import argparse
import ast
import shapefile as sf


def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    ds.close()
    return(XC, YC, Depth)

def read_CTD_points_from_shp(config_dir, model_name):
    shp_file = os.path.join(config_dir,'L2',model_name,'input','dv_shp','CTD')
    r = sf.Reader(shp_file)
    records = r.records()
    shapes = r.shapes()
    ctd_points = []
    for r in range(len(records)):
        ctd_points.append(shapes[r].points[0])
    return(np.array(ctd_points))

def series_to_N_points(series,N):
    #find the total length of the series
    totalDistance=0
    for s in range(len(series[:,0])-1):
        totalDistance+=((series[s,0]-series[s+1,0])**2+(series[s,1]-series[s+1,1])**2)**0.5
    intervalDistance=totalDistance/(N-1)

    #make the list of points
    newSeries=series[0,:]
    currentS = 0
    currentPoint1=series[currentS,:]
    currentPoint2=series[currentS+1,:]
    for p in range(N-2):
        distanceAccrued = 0
        while distanceAccrued<intervalDistance:
            currentLineDistance=((currentPoint1[0]-currentPoint2[0])**2+(currentPoint1[1]-currentPoint2[1])**2)**0.5
            if currentLineDistance<intervalDistance-distanceAccrued:
                distanceAccrued+=currentLineDistance
                currentS+=1
                currentPoint1 = series[currentS, :]
                currentPoint2 = series[currentS + 1, :]
            else:
                distance=intervalDistance-distanceAccrued
                newX=currentPoint1[0]+(distance/currentLineDistance)*(currentPoint2[0]-currentPoint1[0])
                newY = currentPoint1[1] + (distance / currentLineDistance) * (currentPoint2[1] - currentPoint1[1])
                distanceAccrued=intervalDistance+1
                newSeries=np.vstack([newSeries,np.array([newX,newY])])
                currentPoint1=np.array([newX,newY])
    newSeries = np.vstack([newSeries, series[-1,:]])
    return(newSeries)

def read_transect_points_from_shp(config_dir, model_name, transect_name):
    shp_file = os.path.join(config_dir,'L2',model_name,'input','dv_shp',transect_name)
    r = sf.Reader(shp_file)
    shapes = r.shapes()
    points = np.array(shapes[0].points)

    points = series_to_N_points(points,1000)

    XC_boundary = points[:,0]
    YC_boundary = points[:,1]

    return(XC_boundary, YC_boundary)

def create_mask_grid(resolution, L2_XC, L2_YC, L2_Depth, XC_boundary, YC_boundary):

    mask_indices = []
    counter = 1

    mask_grid = np.zeros_like(L2_XC)

    for i in range(len(XC_boundary)):
        if i%int(len(XC_boundary)/100)==0:
            print(i,len(XC_boundary),str(int(100*i/len(XC_boundary)))+'%')
        # dist = ((L2_XC-XC_boundary[i])**2 + (L2_YC-YC_boundary[i])**2)**0.5
        dist = great_circle_distance(XC_boundary[i], YC_boundary[i], L2_XC, L2_YC)
        rows, cols = np.where(dist<resolution*np.sqrt(2))
        for i in range(len(rows)):
            if mask_grid[rows[i],cols[i]]==0 and L2_Depth[rows[i],cols[i]]>0:
                mask_grid[rows[i],cols[i]] = counter
                mask_indices.append([rows[i],cols[i]])
                counter +=1

    mask_indices = np.array(mask_indices)

    return(mask_grid, mask_indices)

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

def create_CTD_mask_grid(L2_XC, L2_YC, L2_Depth, ctd_points):

    mask_indices = []
    counter = 1

    mask_grid = np.zeros_like(L2_XC)

    nonzero_XC = L2_XC[L2_Depth > 0]
    nonzero_YC = L2_YC[L2_Depth > 0]

    for p in range(np.shape(ctd_points)[0]):

        dist_err = 1e22

        dist = great_circle_distance(ctd_points[p,0],ctd_points[p,1],nonzero_XC,nonzero_YC)
        index = np.argmin(dist)

        x = nonzero_XC[index]
        y = nonzero_YC[index]

        ctd_dist = great_circle_distance(ctd_points[p,0], ctd_points[p,1], x, y)

        if ctd_dist<dist_err:

            row,col = np.where(np.logical_and(L2_XC==x,L2_YC==y))
            ctd_row = row[0]
            ctd_col = col[0]
            dist_err = ctd_dist

        print('    - CTD point '+str(p+1)+' ('+str(ctd_points[p,0])+', '+str(ctd_points[p,1])+') will be located at ('+str(ctd_row)+','+
              str(ctd_col)+') with a depth of '+str(L2_Depth[ctd_row,ctd_col])+' (dist err = '+str(dist_err)+')')

        mask_grid[ctd_row,ctd_col] = counter
        mask_indices.append([ctd_row, ctd_col])
        counter += 1

    mask_indices = np.array(mask_indices)

    return(mask_grid, mask_indices)

def output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names_list):
    if output_file_name in os.listdir(output_dir):
        os.remove(os.path.join(output_dir,output_file_name))

    ds = nc4.Dataset(os.path.join(output_dir,output_file_name),'w')

    for m in range(len(mask_names_list)):
        grp = ds.createGroup(mask_names_list[m])
        grp.createDimension('n_points', np.shape(all_mask_dicts[m])[0])
        var = grp.createVariable('source_rows', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 0].astype(int)
        var = grp.createVariable('source_cols', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 1].astype(int)

    ds.close()


########################################################################################################################

def create_dv_masks(config_dir, L2_model_name, print_level):

    # if 'input' not in os.listdir('..'):
    #     os.mkdir(os.path.join('..','input'))
    if 'dv' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'input','dv'))

    ###############################################################################################
    # Read in the grids

    if print_level>=1:
        print('    - Reading in the L2 geometry')

    # read the extended subset of the XC and YC grids
    L2_XC, L2_YC, L2_Depth = read_grid_geometry_from_nc(config_dir, L2_model_name)
    if print_level >= 2:
        print('        - The L2 domain has shape ' + str(np.shape(L2_XC)))

    ###############################################################################################
    # Create the masks

    resolution = 500

    mask_names_list = []
    all_mask_dicts = []

    ctd_points = read_CTD_points_from_shp(config_dir, L2_model_name)

    for boundary in ['Ilulissat_Fjord','Torsukataq_Fjord','CTD']:#,'DJG_Fjord','Sund_Entrance'] #CTD
        if print_level >= 1:
            print('    - Creating the ' + boundary +' mask')

        if boundary!='CTD':
            XC_boundary, YC_boundary = read_transect_points_from_shp(config_dir, L2_model_name, boundary)
            mask_grid, mask_indices = create_mask_grid(resolution, L2_XC, L2_YC, L2_Depth, XC_boundary, YC_boundary)
        else:
            mask_grid, mask_indices = create_CTD_mask_grid(L2_XC, L2_YC, L2_Depth, ctd_points)
        all_mask_dicts.append(mask_indices)
        mask_names_list.append(boundary)

        if print_level>=2:
            print('        - The '+boundary+' mask has '+str(np.shape(mask_indices)[0])+' points')

        # C = plt.imshow(mask_grid,origin='lower')
        # plt.colorbar(C)
        # plt.title(boundary+' ('+str(np.shape(mask_indices)[0])+' points)')
        # plt.show()

        output_file = os.path.join(config_dir, 'L2', L2_model_name, 'input', 'dv', boundary+'_mask.bin')
        mask_grid=np.array(mask_grid)
        mask_grid.ravel(order='C').astype('>f4').tofile(output_file)

    output_dir = os.path.join(config_dir,'L2',L2_model_name,'input')
    output_file_name = 'L2_dv_mask_reference_dict.nc'
    output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names_list)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_dv_masks(config_dir)
