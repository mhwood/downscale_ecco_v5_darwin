
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from pyproj import Transformer
import matplotlib.path as mplPath
import sys

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon


def read_in_domain_boundary(config_dir, model_name):

    grid_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(grid_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    ds.close()

    bottom = np.column_stack([XC[0, :], YC[0, :]])
    right = np.column_stack([XC[:, -1], YC[:, -1]])
    top = np.column_stack([XC[-1, :], YC[-1, :]])
    left = np.column_stack([XC[:,0], YC[:,0]])

    boundary = np.vstack([bottom,right,np.flipud(top),np.flipud(left)])

    boundary = reproject_polygon(boundary, 4326, 3413)

    # plt.plot(boundary[:,0],boundary[:,1])
    # plt.show()

    return(XC, YC, Depth, boundary)


def find_Mankoff_outlets_in_domain(mankoff_dir, domain_boundary):

    file_name = os.path.join(mankoff_dir,'outlets','freshwater','ice','outlets.csv')
    # f = open(file_name)
    # lines = f.read()
    # f.close()
    # lines = lines.split('\n')
    # lines.pop(0)
    outlet_locations = np.genfromtxt(file_name,skip_header=1,delimiter=',')

    points = outlet_locations[:,-2:]
    p = mplPath.Path(domain_boundary)
    inside = p.contains_points(points)

    domain_outlet_locations = outlet_locations[inside,:]
    domain_outlet_locations = domain_outlet_locations[domain_outlet_locations[:,5]<=-10,:]

    # plt.plot(domain_boundary[:, 0], domain_boundary[:, 1])
    # plt.plot(domain_outlet_locations[:,-2],domain_outlet_locations[:,-1],'k.')
    # plt.show()

    outlet_points_xyz = np.column_stack([domain_outlet_locations[:,-4:-2],domain_outlet_locations[:,5]])
    outlet_points_coast_id = domain_outlet_locations[:,0]

    return(outlet_points_xyz, outlet_points_coast_id)


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


def map_Mankoff_outlets_to_model_domain(outlet_points_xyz, XC_faces, YC_faces, Depth_faces, hFac_faces, RC, print_level):

    printing = True

    model_fijxyz = np.zeros((np.shape(outlet_points_xyz)[0],6))

    for i in range(np.shape(outlet_points_xyz)[0]):

        if print_level>=5:
            print('    - Point ' + str(i))
            print('          - Discharge Location : ' + str(outlet_points_xyz[i, 0]) + ', ' + str(outlet_points_xyz[i, 1]) + ' at depth ' + str(outlet_points_xyz[i, 2]) + ' m')

        point_err=1e10

        for face in list(XC_faces.keys()):
            domain_X = XC_faces[face]
            domain_Y = YC_faces[face]
            domain_depth = Depth_faces[face]

            xyz = np.column_stack([domain_X.ravel(), domain_Y.ravel(), domain_depth.ravel()])

            subset_xyz = xyz[xyz[:, 2] >= 0.5 * outlet_points_xyz[i, 2], :]
            # subset_xyz = sub

            distance = great_circle_distance(lon_ref = outlet_points_xyz[i,0], lat_ref = outlet_points_xyz[i,1],
                                             Lon = subset_xyz[:,0], Lat = subset_xyz[:,1])
            index = np.argmin(distance)

            if distance[index]<point_err:
                point_err = distance[index]

                row, col = np.where(np.logical_and(domain_X==subset_xyz[index,0],domain_Y==subset_xyz[index,1]))

                k = np.argmin(np.abs(RC-outlet_points_xyz[i,2]))
                #  this ensures that the discharge will be put in a cell that's atleast 75% wet and at least 200 m deep
                for n in range(np.shape(hFac_faces[face])[0]):
                    if hFac_faces[face][k,row[0],col[0]]<0.75 and k!=0:
                        k-=1

                model_fijxyz[i, 0] = face
                model_fijxyz[i, 1] = row[0]
                model_fijxyz[i, 2] = col[0]

                model_fijxyz[i, 3] = subset_xyz[index,0]
                model_fijxyz[i, 4] = subset_xyz[index,1]
                model_fijxyz[i, 5] = RC[k]#subset_xyz[index,2]

        if print_level >= 5:
            print('          - Closest point thats deep enough: ' + str(subset_xyz[index, 0]) + ', ' + str(subset_xyz[index, 1]) + ', at depth ' + str(subset_xyz[index, 2]))
            print('                at a distance of ' + str(distance[index]) + ' m')

    # # plot some lines to see how far away the source and closest model cell are
    # for i in range(np.shape(model_fijxyz)[0]):
    #     print(model_fijxyz[i,:])
    #     plt.plot([outlet_points_xyz[i,0],model_ijxyzd[i,2]],[outlet_points_xyz[i,1],model_ijxyzd[i,3]],'k-')
    # plt.show()

    # plt.pcolormesh(domain_X,domain_Y,domain_depth,alpha=0.5)
    # plt.plot(model_ijxyzd[:,2],model_ijxyzd[:,3],'k.')
    # plt.plot(outlet_points_xyz[:,0],outlet_points_xyz[:,1],'g.')
    # plt.show()

    return(model_fijxyz)


def create_iceplume_mask(Lf, config_dir, model_name, outlet_points_xyz, model_fijxyz, hFac_faces, sNx, sNy):

    mask_faces = {}
    for face in list(hFac_faces.keys()):
        mask_faces[face] = np.zeros_like(hFac_faces[face])

    points_with_multiple_sources = 0
    for i in range(np.shape(model_fijxyz)[0]):
        face = model_fijxyz[i,0].astype(int)
        row = model_fijxyz[i,1].astype(int)
        col = model_fijxyz[i,2].astype(int)
        k = model_fijxyz[i,3].astype(int)
        if mask_faces[face][k, row, col]==0:
            mask_faces[face][k, row, col] = 3
        else:
            points_with_multiple_sources += 1

    print('        - '+str(points_with_multiple_sources)+' will receive fluxes from multiple sources ('+str(np.shape(model_fijxyz)[0])+' total sources)')

    mask_compact = Lf.read_faces_to_compact(mask_faces, sNx, sNy)

    if 'iceplume' not in os.listdir(os.path.join(config_dir, 'L3', model_name,'input')):
        os.mkdir(os.path.join(config_dir, 'L3', model_name,'input','iceplume'))

    output_file = os.path.join(config_dir, 'L3', model_name, 'input', 'iceplume', '_'.join(model_name.split('_')[:2])+'_iceplume_mask.bin')
    mask_compact.ravel(order='C').astype('>f4').tofile(output_file)

    output_lines = 'N,Mankoff_X,Mankoff_Y,Mankoff_Depth,Model_X,Model_Y,Model_Depth,Model_Face,Model_Row,Model_Col'
    for i in range(np.shape(model_fijxyz)[0]):
        output_lines += '\n'+str(i+1)+','
        output_lines += str(outlet_points_xyz[i, 0]) + ',' + str(outlet_points_xyz[i, 1]) + ',' + str(outlet_points_xyz[i, 2])+','
        output_lines += str(model_fijxyz[i, 3]) + ',' + str(model_fijxyz[i, 4]) + ',' + str(model_fijxyz[i, 5]) + ','
        output_lines += str(model_fijxyz[i, 0]) + ',' + str(model_fijxyz[i, 1]) + ',' + str(model_fijxyz[i, 2])

    f = open(os.path.join(config_dir, 'L3', model_name, 'input', 'iceplume', '_'.join(model_name.split('_')[:2])+'_iceplume_mask_ref.csv'),'w')
    f.write(output_lines)
    f.close()



def create_annual_iceplume_files(config_dir, model_name, mankoff_dir, domain_X, domain_Y, model_fijxyz, outlet_points_coast_id, years):


    for year in years:
        print('        - Creating the file for year '+str(year))

        # find out how many days will be in the output file
        if year%4==0:
            n_days = 366
        else:
            n_days = 365

        # make output grid faces for this file
        output_grid = np.zeros_like(domain_Y)

        # read in the Mankoff data for this year
        discharge_file = os.path.join(mankoff_dir,'runoff','ice','runoff','MAR_'+str(year)+'.nc')
        ds = nc4.Dataset(discharge_file)
        coast_id = ds.variables['station'][:]
        runoff = ds.variables['runoff'][:,:]
        ds.close()

        for i in range(len(outlet_points_coast_id)):
            # face = model_fijxyz[i, 0].astype(int)
            # row = model_fijxyz[i, 1].astype(int)
            # col = model_fijxyz[i, 2].astype(int)


            id = int(outlet_points_coast_id[i])
            id_index = np.where(coast_id==id)[0][0]

            plt.plot(runoff[id_index,:])
            plt.title(str(i)+','+str(coast_id))
            plt.show()

            # output_faces[face][:,row,col] += runoff[id_index,:].ravel()

        # # C = plt.pcolormesh(output_grid)
        # # plt.colorbar(C)
        # # plt.show()

        # output_compact = Lf.read_faces_to_compact(output_faces, sNx, sNy)
        #
        # output_file = os.path.join(config_dir, 'L3', model_name, 'input', 'iceplume', '_'.join(model_name.split('_')[:2])+'_Qsg_'+str(year))
        # output_compact.ravel(order='C').astype('>f4').tofile(output_file)



def create_L3_iceplume_files(config_dir, model_name, mankoff_dir, years, print_level):

    if print_level >= 1:
        print('    - Creating the iceplume files for the '+model_name+' model from Mankoff subglacial discharge data')

    if print_level >= 1:
        print('    - Reading in the model boundary')
    domain_X, domain_Y, domain_depth, domain_boundary = read_in_domain_boundary(config_dir, model_name)

    # step 1: find all of the Mankoff outlets in the domain
    if print_level >= 1:
        print('    - Identifying subglacial discharge outlets within the domain from Mankoff references')
    outlet_points_xyz, outlet_points_coast_id = find_Mankoff_outlets_in_domain(mankoff_dir, domain_boundary)


    # plt.plot(domain_boundary[:,0], domain_boundary[:,1],'k-')
    # plt.plot(outlet_points_xyz[:,0],outlet_points_xyz[:,1],'g.')
    # plt.show()

    # # make sure all the depths are positive
    # RC = np.abs(RC)
    # outlet_points_xyz[:,2] = np.abs(outlet_points_xyz[:,2])

    # # step 2: create a 2D mask for the domain
    # if print_level >= 1:
    #     print('    - Mapping outlet locations to model grid cells')
    # model_fijxyz = map_Mankoff_outlets_to_model_domain(outlet_points_xyz, XC_faces, YC_faces, Depth_faces, hFac_faces, RC, print_level)
    model_fijxyz = []
    #
    # # step 3: create iceplume mask
    # if print_level >= 1:
    #     print('    - Creating a mask for use in the iceplume package')
    # create_iceplume_mask(Lf, config_dir, model_name, outlet_points_xyz, model_fijxyz, hFac_faces, sNx, sNy)

    # step 4: create the iceplume file
    if print_level >= 1:
        print('    - Creating a file with a timeseries of fluxes for use in the iceplume package')
    create_annual_iceplume_files(config_dir, model_name, mankoff_dir, domain_X, domain_Y, model_fijxyz, outlet_points_coast_id, years)




