
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from pyproj import Proj, Transformer
import matplotlib.path as mplPath
# from osgeo import gdal
import shapefile
import datetime
import argparse

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1):
    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))
    # inProj = Proj(init='epsg:'+str(inputCRS))
    # outProj = Proj(init='epsg:'+str(outputCRS))
    x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon


def get_model_grid_boundary(config_dir,model_name):

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',model_name+'_grid.nc'))
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    rA = ds.variables['rA'][:, :]
    # rA = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    ds.close()

    left = np.column_stack([XC[:, 0], YC[:, 0]])
    right = np.column_stack([XC[:, -1], YC[:, -1]])
    bottom = np.column_stack([XC[0, :], YC[0, :]])
    top = np.column_stack([XC[-1, :], YC[-1, :]])

    boundary = np.vstack([bottom,right,np.flipud(top),np.flipud(left)])

    boundary_3413 = reproject_polygon(boundary,4326,3413)#,x_column=1,y_column=0)

    points = np.column_stack([XC.ravel(), YC.ravel()])
    points = reproject_polygon(points,4326,3413)
    X = np.reshape(points[:,0],np.shape(XC))
    Y = np.reshape(points[:,1],np.shape(YC))

    return(X, Y, rA, Depth, boundary, boundary_3413)


def get_land_outlet_locations(runoff_dir, model_boundary_3413):

    file_path = os.path.join(runoff_dir,'Mankoff_liquid','freshwater','land','outlets.csv')

    f = open(file_path)
    lines = f.read()
    f.close()
    header = lines.split('\n')[0].split(',')

    points = np.genfromtxt(file_path,delimiter=',',skip_header=1)

    p = mplPath.Path(model_boundary_3413)
    inside = p.contains_points(points[:,3:5])
    points = points[inside,:]

    return(header,points)


def get_ice_outlet_locations(runoff_dir, model_boundary_3413):

    file_path = os.path.join(runoff_dir,'Mankoff_liquid','freshwater','ice','outlets.csv')

    f = open(file_path)
    lines = f.read()
    f.close()
    header = lines.split('\n')[0].split(',')

    points = np.genfromtxt(file_path,delimiter=',',skip_header=1)

    p = mplPath.Path(model_boundary_3413)
    inside = p.contains_points(points[:,9:11])
    points = points[inside,:]

    return(header,points)


def get_ice_fronts_in_domain(termpicks_file, glacier_IDs):

    fronts = []
    front_IDs = []
    sf = shapefile.Reader(termpicks_file)
    records = sf.records()
    shapes = sf.shapes()
    for r in range(len(records)):
        record = records[r]
        if int(record[0]) in glacier_IDs:
            if int(record[2])>2000:
                fronts.append(np.array(shapes[r].points))
                front_IDs.append(record[0])

    return(fronts, front_IDs)


def categorize_points(points,fronts,front_IDs,x_col,y_col):

    center_points = np.zeros((len(fronts),2))
    for f in range(len(fronts)):
        center_points[f, 0] = np.mean(fronts[f][:, 0])
        center_points[f, 1] = np.mean(fronts[f][:, 1])

    threshold = 1000

    categories = np.zeros((np.shape(points)[0],))

    for p in range(np.shape(points)[0]):
        point_dist = (center_points[:,0]-points[p,x_col])**2 + (center_points[:,1]-points[p,y_col])**2
        front_index = np.argmin(point_dist)
        front = fronts[front_index]
        front_dist = ((front[:,0]-points[p,x_col])**2 + (front[:,1]-points[p,y_col])**2)**0.5
        if np.min(front_dist)<threshold:
            categories[p] = front_IDs[front_index]

    return(categories)


def read_discharge_from_file(runoff_dir, subset, discharge_type):

    print('                - Reading in the big file...')
    ds = nc4.Dataset(os.path.join(runoff_dir,'Mankoff_liquid','freshwater',discharge_type,'MAR.nc'))
    time = ds.variables['time'][:]
    station = ds.variables['station'][:]
    discharge = ds.variables['discharge'][:,:]
    ds.close()

    print('                - Looping through the points')
    glacier_discharge = np.zeros((len(time),))
    for i in range(np.shape(subset)[0]):
        cat = int(subset[i, 0])
        index = np.where(station == cat)[0]
        if len(index)>0:
            index = index[0]
            index_subset = discharge[index,:]
            glacier_discharge[~np.isnan(index_subset)] += index_subset[~np.isnan(index_subset)]
        else:
            print('                - Missing '+discharge_type+' point with cat = '+str(cat))

    del discharge

    return(time,glacier_discharge)


def write_discharge_timeseries_to_nc(project_dir,glacierID,time,ice_discharge,land_discharge):

    ds = nc4.Dataset(os.path.join(project_dir,'Glacier_'+str(glacierID)+'_discharge.nc'), 'w')
    ds.createDimension('time',len(time))

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = time
    tvar.note = 'days since 1950/1/1'

    ivar = ds.createVariable('ice_discharge', 'f4', ('time',))
    ivar[:] = ice_discharge
    ivar.note = 'liquid discharge that first came from ice'

    lvar = ds.createVariable('land_discharge', 'f4', ('time',))
    lvar[:] = land_discharge
    lvar.note = 'liquid discharge that first came from land'

    ds.close()


def generate_glacier_subglacial_discharge_files(project_dir, runoff_dir, glacierIDs, land_points, land_categories, ice_points, ice_categories):

    if 'glaciers' not in os.listdir(project_dir):
        os.mkdir(os.path.join(project_dir,'glaciers'))

    for glacierID in glacierIDs:

        if 'Glacier_'+str(glacierID)+'_discharge.nc' not in \
                 os.listdir(os.path.join(project_dir,'glaciers')):

            print('        - Working on glacier '+str(glacierID))
            land_subset = land_points[land_categories==glacierID,:]
            ice_subset = ice_points[ice_categories==glacierID,:]

            print('            - Glacier '+str(glacierID)+' has '+str(np.shape(land_subset)[0])+
                  ' land points and '+str(np.shape(ice_subset)[0])+' ice points')

            print('            - Working on the ice discharge timeseries...')
            time, ice_discharge = read_discharge_from_file(runoff_dir, ice_subset, discharge_type='ice')

            print('            - Working on the land discharge timeseries...')
            _, land_discharge = read_discharge_from_file(runoff_dir, land_subset, discharge_type='land')

            print('            - Writing out to nc...')
            write_discharge_timeseries_to_nc(project_dir, glacierID, time, ice_discharge, land_discharge)


def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)


def time_to_dec_yr(time):

    dec_yr = np.zeros((len(time),))
    for t in range(len(time)):
        date = datetime.datetime(1950,1,1)+datetime.timedelta(days=int(time[t]))
        dec_yr[t] = YMD_to_DecYr(date.year, date.month, date.day)

    return(dec_yr)


def map_discharge_points_to_model_points(X, Y, Depth, points,x_col,y_col):

    model_points = np.column_stack([X.ravel(),Y.ravel()])
    model_points = model_points[Depth.ravel()!=0,:]

    rows = np.zeros((np.shape(points)[0],))
    cols = np.zeros((np.shape(points)[0],))

    for p in range(len(rows)):
        dist = ((model_points[:,0]-points[p,x_col])**2 + (model_points[:,1]-points[p,y_col])**2)**0.5
        index = np.where(dist==np.min(dist))[0][0]
        Dist = (X-model_points[index,0])**2 + (Y-model_points[index,1])**2
        r, c = np.where(Dist==np.min(Dist))
        # print('        point '+str(points[p,x_col])+', '+str(points[p,y_col])+
        #       ' maps to '+str(X[r[0],c[0]])+', '+str(Y[r[0],c[0]])+' (dist = '+str(np.min(dist))+')')
        rows[p] = r[0]
        cols[p] = c[0]

    return(rows,cols)


def read_runoff_from_file(runoff_dir, discharge_type, year_indices, X, Y,
                          points, categories, model_rows, model_cols):

    runoff = np.zeros((np.sum(year_indices), np.shape(X)[0], np.shape(X)[1]))

    print('                - Reading in the big file...')
    ds = nc4.Dataset(os.path.join(runoff_dir,'Mankoff_liquid','freshwater',discharge_type, 'MAR.nc'))
    time = ds.variables['time'][:]
    station = ds.variables['station'][:]
    discharge = ds.variables['discharge'][:,:]
    ds.close()

    print('                - Looping through the points')
    for i in range(np.shape(points)[0]):
        if categories[i]==0:
            cat = int(points[i, 0])
            index = np.where(station == cat)[0]
            if len(index)>0:
                index = index[0]
                index_subset = discharge[index,:]
                index_subset = index_subset[year_indices]
                index_subset[np.isnan(index_subset)] = 0
                runoff[:,int(model_rows[i]),int(model_cols[i])] += index_subset
            else:
                print('        - Missing '+discharge_type+' point with cat = '+str(cat))

    del discharge

    return(runoff)

def generate_annual_runoff_files(runoff_dir, config_dir,model_name,  years, X, Y, cell_area,
                                 land_points, land_categories, land_rows, land_cols,
                                 ice_points, ice_categories, ice_rows, ice_cols):

    ds = nc4.Dataset(os.path.join(runoff_dir,'Mankoff_liquid','freshwater','land', 'MAR.nc'))
    time = ds.variables['time'][:]
    ds.close()

    for year in years:
        dec_yr = time_to_dec_yr(time)
        year_indices = np.logical_and(dec_yr >= year, dec_yr < year + 1)

        land_runoff = read_runoff_from_file(runoff_dir, 'land', year_indices, X, Y,
                                            land_points, land_categories, land_rows, land_cols)

        print(np.shape(land_runoff))

        ice_runoff = read_runoff_from_file(runoff_dir, 'ice', year_indices, X, Y,
                                           ice_points, ice_categories, ice_rows, ice_cols)

        total_runoff = land_runoff + ice_runoff

        # divide by the area of each cell to provide runoff in m/s
        for day in range(np.shape(total_runoff)[0]):
            day_runoff = total_runoff[day, :, :]/cell_area
            total_runoff[day, :, :] = day_runoff

        print('            - Output stats: min: ' + str(np.min(total_runoff)) + ', max: ' + str(np.max(total_runoff)))

        output_file = os.path.join(config_dir,'L2',model_name,'input','exf','L2_exf_RUNOFF_Mankoff_'+str(year))
        total_runoff.ravel(order='C').astype('>f4').tofile(output_file)



def create_L2_iceplume_files(config_dir, model_name, runoff_dir, termpicks_file, glacier_IDs, years, print_level):

    if print_level>=1:
        print('    - Creating the iceplume files for the '+model_name+' model')

    print('    - Reading model geometry')
    X, Y, cell_area, Depth, model_boundary, model_boundary_3413 = get_model_grid_boundary(config_dir,model_name)

    print('    - Reading land and ice points')
    land_column_names, land_points = get_land_outlet_locations(runoff_dir, model_boundary_3413)
    ice_column_names, ice_points = get_ice_outlet_locations(runoff_dir, model_boundary_3413)

    # output_path = os.path.join(project_dir,'Figures','Ocean','Discharge Strategy','Mankoff Land Points.png')
    # create_land_points_plot(output_path,project_dir,land_points)
    #
    # output_path = os.path.join(project_dir,'Figures','Ocean','Discharge Strategy','Mankoff Ice Points.png')
    # create_ice_points_plot(output_path,project_dir,ice_points)

    print('    - Getting fronts')
    fronts, front_IDs = get_ice_fronts_in_domain(termpicks_file, glacier_IDs)

    # output_path = os.path.join(project_dir,'Figures','Ocean','Discharge Strategy','TermPicks Ice Fronts.png')
    # create_ice_fronts_plot(output_path, project_dir, fronts)

    print('    - Categorizing the land points')
    land_categories = categorize_points(land_points,fronts,front_IDs,x_col=3,y_col=4)

    print('    - Categorizing the ice points')
    ice_categories = categorize_points(ice_points,fronts,front_IDs,x_col=9,y_col=10)

    # output_path = os.path.join(project_dir,'Figures','Ocean','Discharge Strategy','Mankoff Land Point Categories.png')
    # create_land_points_categories_plot(output_path, project_dir, land_points, land_categories)
    #
    # output_path = os.path.join(project_dir,'Figures','Ocean','Discharge Strategy','Mankoff Ice Point Categories.png')
    # create_ice_points_categories_plot(output_path, project_dir, ice_points, ice_categories)

    # print('    - Generating the subglacial discharge files')
    # output_dir = os.path.join(config_dir, model_level, model_name, 'input')
    # if 'glacier_sgd' not in os.listdir(output_dir):
    #     os.mkdir(os.path.join(output_dir, 'glacier_sgd'))
    # output_dir = os.path.join(output_dir, 'glacier_sgd')
    # generate_glacier_subglacial_discharge_files(output_dir, runoff_dir, glacier_IDs, land_points, land_categories, ice_points, ice_categories)

    # print('    - Plotting an example discharge timeseries')
    # output_path = os.path.join(project_dir,'Figures','Ocean','Discharge Strategy','Daugaard Jensen Discharge.png')
    # plot_subglacial_discharge_for_glacier(output_path, project_dir,glacierID=118,year=2000)

    print('    - Mapping runoff points to model points')
    land_rows, land_cols = map_discharge_points_to_model_points(X, Y, Depth, land_points, x_col=3, y_col=4)
    ice_rows, ice_cols = map_discharge_points_to_model_points(X, Y, Depth, ice_points, x_col=9, y_col=10)

    print('    - Generating the runoff files')
    generate_annual_runoff_files(runoff_dir, config_dir, model_name, years, X, Y, cell_area,
                                 land_points, land_categories, land_rows, land_cols,
                                 ice_points, ice_categories, ice_rows, ice_cols)

    # print('    - Plotting an example runoff field')
    # output_path = os.path.join(project_dir,'Figures','Ocean','Discharge Strategy','Runoff Field Example.png')
    # plot_runoff_field(output_path, config_dir, X, Y, year=2000, year_day=183)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-l", "--model_level", action="store",
                        help="The level of the model (e.g. L2).", dest="model_level",
                        type=str, required=True)

    parser.add_argument("-m", "--model_name", action="store",
                        help="The name of the model (e.g. L2_Disko_Bay).", dest="model_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    model_name = args.model_name
    model_level = args.model_level

    create_L2_iceplume_files(config_dir, model_level, model_name, print_level=5)
