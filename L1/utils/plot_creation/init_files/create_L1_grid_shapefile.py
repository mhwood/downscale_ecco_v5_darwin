
import os
import numpy as np
import netCDF4 as nc4
import shapefile
import matplotlib.pyplot as plt
import argparse


def grid_to_shp(config_dir,model_name):
    # grid_file = '/Volumes/mhwood/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
    #             'downscaled_greenland/nc_grids/L3_Scoresby_Sund_grid.nc'
    grid_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')

    ds = nc4.Dataset(grid_file)
    XG = ds.variables['XG'][:,:]
    YG = ds.variables['YG'][:,:]
    # Depth = ds.variables['Depth'][:,:]
    ds.close()

    if model_name=='L1_mac_delta':
        XG-=360

    output_file = os.path.join(config_dir,'shp_grids',model_name+'_grid')

    sf = shapefile.Writer(output_file)

    sf.field('row','N')
    sf.field('col','N')
    # sf.field('depth','N')

    for row in range(np.shape(XG)[0]-1):
        for col in range(np.shape(XG)[1]-1):
            # if int(Depth[row,col])>0:
            # polygon = [[XG[row,col]-XG[row,col],YG[row,col]-YG[row,col]],
            #            [XG[row,col]+XG[row,col],YG[row,col]-YG[row,col]],
            #            [XG[row,col]+XG[row,col],YG[row,col]+YG[row,col]],
            #            [XG[row,col]-XG[row,col],YG[row,col]+YG[row,col]]]
            polygon = [[XG[row, col], YG[row, col]],
                       [XG[row+1, col], YG[row+1, col]],
                       [XG[row+1, col+1], YG[row+1, col+1]],
                       [XG[row, col+1], YG[row, col+1]]]
            sf.poly([polygon])
            sf.record(row,col)#,int(Depth[row,col]))


    sf.close()

    f = open(output_file + '.prj', 'w')
    f.write(
        'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]')
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-m", "--model_name", action="store",
                        help="The name of the model (e.g. L1_W_Greenland.", dest="model_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    model_name = args.model_name

    grid_to_shp(config_dir,model_name)