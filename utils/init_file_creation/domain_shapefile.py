
import os
import numpy as np
import netCDF4 as nc4
import shapefile
import matplotlib.pyplot as plt
import argparse

def create_shp(config_dir, model_level, model_name, print_level):

    if print_level>=1:
        print('    - Creating the shapefile for the '+model_name+' model')

    # grid_file = '/Volumes/mhwood/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
    #             'downscaled_greenland/nc_grids/L3_Scoresby_Sund_grid.nc'
    grid_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')

    ds = nc4.Dataset(grid_file)
    XG = ds.variables['XG'][:,:]
    YG = ds.variables['YG'][:,:]
    Depth = ds.variables['Depth'][:,:]
    ds.close()

    # output_file = '/Users/michwood/Documents/Research/Projects/Scoresby Sund/Map/Shapefiles/Domains/L3_Scoresby_Sund'
    output_dir = os.path.join(config_dir,model_level,model_name,'input')
    if 'domain_shp' not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,'domain_shp'))
    output_dir = os.path.join(output_dir,'domain_shp')
    output_file = os.path.join(output_dir,model_name)

    sf = shapefile.Writer(output_file)

    sf.field('row','N')
    sf.field('col','N')
    sf.field('depth','N')

    for row in range(np.shape(XG)[0]-1):
        for col in range(np.shape(XG)[1]-1):
            if int(Depth[row,col])>0:
                # polygon = [[XG[row,col]-XG[row,col],YG[row,col]-YG[row,col]],
                #            [XG[row,col]+XG[row,col],YG[row,col]-YG[row,col]],
                #            [XG[row,col]+XG[row,col],YG[row,col]+YG[row,col]],
                #            [XG[row,col]-XG[row,col],YG[row,col]+YG[row,col]]]
                polygon = [[XG[row, col], YG[row, col]],
                           [XG[row+1, col], YG[row+1, col]],
                           [XG[row+1, col+1], YG[row+1, col+1]],
                           [XG[row, col+1], YG[row, col+1]]]
                sf.poly([polygon])
                sf.record(row,col,int(Depth[row,col]))

    sf.close()

    prj_file = output_file+'.prj'
    f = open(prj_file,'w')
    f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],'
            'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
    f.close()

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

    create_shp(config_dir, model_level, model_name, print_level=5)
