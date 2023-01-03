
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
import shutil
import cmocean.cm as cm
import argparse
from matplotlib.patches import Rectangle
import datetime as dt
from datetime import datetime, timedelta
from osgeo import gdal
import sys


def read_field_from_monthly_ncs(config_dir, config_name, results_dir, field_name):

    if field_name in ['SIarea','SIheff','SIhsnow','SIuice','SIvice']:
        file_prefix = 'SI_daily_snap'
    elif field_name in ['Theta','Salt']:
        file_prefix = 'TS_surf_daily_snap'
    elif field_name in ['Vorticity','Speed','Uvel','Vvel']:
        file_prefix = 'vel_surf_daily_snap'
    elif field_name in ['Vorticity_AW']:
        file_prefix = 'vel_AW_daily_snap'
    elif field_name in ['Theta_AW','Salt_AW']:
        file_prefix = 'TS_AW_daily_snap'
    elif field_name in ['EXFaqh','EXFatemp','EXFpreci','EXFroff','EXFlwdn','EXFswdn','EXFuwind','EXFvwind']:
        file_prefix = 'EXF_day_snap'
    elif field_name in ['DIC','FeT','NH4','NO2','NO3','PO4','SiO2']:
        file_prefix = 'BGC_daily_consts'
    elif field_name in ['Chl01','Chl02','Chl03','Chl04','Chl05']:
        file_prefix = 'BGC_daily_Chl'
    elif field_name in ['c01','c02','c03','c04','c05','c06','c07']:
        file_prefix = 'BGC_daily_cx'
    elif field_name in ['DOC','DOFe','DON','DOP']:
        file_prefix = 'BGC_daily_DO'
    elif field_name in ['Alk','O2','PIC']:
        file_prefix = 'BGC_daily_misc'
    elif field_name in ['POC','POFe','PON','POP','POSi']:
        file_prefix = 'BGC_daily_PO'
    elif field_name in ['EtaN']:
        file_prefix = 'EtaN_day_snap'
    else:
        raise ValueError('Variable not yet implemented')

    grid_started = False

    year_months = []
    for file_name in os.listdir(os.path.join(results_dir, file_prefix)):
        if file_name[-3:]=='.nc' and file_name[0]!='.':
            year_month = file_name.split('.')[1]
            if year_month not in year_months:
                year_months.append(year_month)

    if field_name=='Vorticity' or field_name=='Vorticity_AW':
        ds = nc4.Dataset(os.path.join(config_dir, 'nc_grids', config_name + '_grid.nc'))
        YC = ds.variables['YC'][:, :]
        dxC = ds.variables['dxC'][:, :]
        dyC = ds.variables['dyC'][:, :]
        rAz = ds.variables['rAz'][:, :]
        ds.close()

        dxC = dxC[:, 1:]
        dyC = dyC[1:, :]
        rAz = rAz[1:, 1:]

        Omega = 7.2921e-5
        f_grid = 2 * Omega * np.sin(np.deg2rad(YC))

    # numerator = np.diff(uvel * L1_DXC, axis=0)[:, :-1] + np.diff(vvel * L1_DYC, axis=1)[:-1, :]
    # denominator = L1_RAZ[:-1, :-1]
    # zeta = np.zeros_like(numerator)
    # zeta[denominator != 0] = numerator[denominator != 0] / denominator[denominator != 0]
    # var_field = zeta / L1_f_grid[:-1, :-1]

    # year_months = ['199201']#,'199202','199203','199204','199205','199206',
    #                #'199207','199208','199209','199210','199211','199212']

    for year_month in year_months:
        for file_name in os.listdir(os.path.join(results_dir,file_prefix)):
            if file_name.split('.')[1]==year_month:
                print('    - Reading from '+file_name)

                if field_name in ['Speed']:
                    ds = nc4.Dataset(os.path.join(results_dir, file_prefix, file_name))
                    uvel = ds.variables['Uvel'][:, :, :]
                    vvel = ds.variables['Vvel'][:, :, :]
                    iter_numbers = ds.variables['iterations'][:]
                    ds.close()
                    field = (uvel ** 2 + vvel ** 2) ** 0.5
                elif field_name in ['Vorticity','Vorticity_AW']:
                    ds = nc4.Dataset(os.path.join(results_dir, file_prefix, file_name))
                    uvel = ds.variables['Uvel'][:, :, :]
                    vvel = ds.variables['Vvel'][:, :, :]
                    iter_numbers = ds.variables['iterations'][:]
                    ds.close()
                    field = np.zeros((np.shape(uvel)[0],np.shape(uvel)[1]-1,np.shape(uvel)[2]-1))
                    for time_step in range(np.shape(uvel)[0]):
                        numerator = np.diff(uvel[time_step,:,:] * dxC, axis=0)[:, :-1] + np.diff(vvel[time_step,:,:] * dyC, axis=1)[:-1, :]
                        denominator = rAz[:-1, :-1]
                        zeta = np.zeros_like(numerator)
                        zeta[denominator != 0] = numerator[denominator != 0] / denominator[denominator != 0]
                        field[time_step,:,:] = zeta / f_grid[:-1, :-1]
                else:
                    ds = nc4.Dataset(os.path.join(results_dir, file_prefix, file_name))
                    field = ds.variables[field_name.split('_')[0]][:, :, :]
                    iter_numbers = ds.variables['iterations'][:]
                    ds.close()
                    field = field[:, : ,:]

                if not grid_started:
                    grid_started = True
                    grid = field
                    all_iter_numbers = np.reshape(iter_numbers,(np.size(iter_numbers),1))
                else:
                    grid = np.concatenate([grid,field],axis=0)
                    all_iter_numbers = np.concatenate([all_iter_numbers,np.reshape(iter_numbers,(np.size(iter_numbers),1))],axis=0)

    return (grid, all_iter_numbers)

def read_background_imagery(file_path):
    ds = gdal.Open(file_path)
    R = np.array(ds.GetRasterBand(1).ReadAsArray())
    G = np.array(ds.GetRasterBand(2).ReadAsArray())
    B = np.array(ds.GetRasterBand(3).ReadAsArray())
    rows = np.shape(R)[0]
    cols = np.shape(R)[1]
    R = R.reshape((rows,cols, 1))
    G = G.reshape((rows, cols, 1))
    B = B.reshape((rows, cols, 1))
    image = np.concatenate([B,G,R],axis=2)
    brightness_factor = 0.1 # 0 - 1
    image = (np.max(image)-np.max(image)*(brightness_factor))/(np.max(image)-np.min(image))*(image-np.min(image))+np.max(image)*(brightness_factor)
    image[image<0]=0
    image[image>1]=1
    transform = ds.GetGeoTransform()
    extents = [transform[0],transform[0]+transform[1]*np.shape(image)[1],transform[3]+ transform[5] * np.shape(image)[0], transform[3]]
    # x_resolution = transform[1]
    # y_resolution = transform[5]
    return(image,extents)

def iter_number_to_date(iter_number):
    seconds_per_iter = 300
    total_seconds = iter_number*seconds_per_iter
    date = datetime(1992,1,1) + timedelta(seconds=total_seconds)
    return(date)

def date_to_iter_number(date):
    seconds_per_iter = 300
    total_seconds = (date-datetime(1992,1,1)).total_seconds()
    iter_number = total_seconds/seconds_per_iter
    return(iter_number)

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = dt.datetime(year,month,day,hour,minute,second)
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def create_panel_plot(output_dir, file_name, metadata_dict, field_name, field_grid, i, iter_number, land_mask, aw_depth_mask, depth, add_background_imagery):

    ############################################################################
    # get and generate some metadat information

    metadata = metadata_dict[field_name]
    vmin = metadata[0]
    vmax = metadata[1]
    cmap = metadata[2]
    units = metadata[3]
    plot_title = metadata[4]

    if np.any(field_grid!=0):
        vmin = np.min(field_grid[field_grid != 0])
        vmax = np.max(field_grid[field_grid != 0])
        print('  vmin: '+str(vmin)+', vmax: '+str(vmax))

    date = iter_number_to_date(iter_number)
    year = date.year
    dec_yr = YMD_to_DecYr(date.year, date.month, date.day)
    min_iter = date_to_iter_number(datetime(date.year, 1, 1))
    max_iter = date_to_iter_number(datetime(date.year + 1, 1, 1))

    if add_background_imagery:
        raise ValueError("Need to generate the background imagery first before using this option")
        file_path = ""
        background_image, extents = read_background_imagery(file_path)

    ############################################################################
    # make the plot

    fig = plt.figure(figsize=(9, 8))
    plt.style.use('dark_background')

    gs2 = GridSpec(17, 12, left=0.05, right=0.95, top = 0.95, bottom = 0.05, hspace=0.05)

    ax1 = fig.add_subplot(gs2[:-2, :-1])

    if add_background_imagery:
        ax1.imshow(background_image, extent=extents, alpha=0.7)

    plot_grid = np.copy(field_grid[i, :, :])
    if 'AW' in field_name:
        plot_grid = np.ma.masked_where(plot_grid==0,plot_grid)
        # land_mask = np.ma.masked_where(land_mask,plot_grid==0)
    else:
        plot_grid = np.ma.masked_where(depth == 0, plot_grid)
    C = ax1.imshow(plot_grid, origin='lower', vmin=vmin, vmax=vmax, cmap=cmap)

    if 'AW' in field_name and add_background_imagery:
        ax1.imshow(aw_depth_mask, origin='lower', vmin=0, vmax=1, cmap='Greys')
        ax1.contour(depth,levels=[0.1],colors='w',linewidths=0.4)

    if not add_background_imagery:
        ax1.imshow(land_mask, origin='lower', vmin=0, vmax=2, cmap='Greys')

    ax1.set_xticks([])
    ax1.set_yticks([])

    cbar = plt.colorbar(C, fraction=0.026, pad=0.04)
    cbar.set_label(units, rotation=90)
    plt.title(plot_title)

    # ax1.text(5, 5, 'Timestep: ' + str(int(iter_number)), ha='left', va='bottom', color='k')

    min_year = 1992
    max_year = 2022
    ax2 = fig.add_subplot(gs2[-1, 1:-2])
    width = (dec_yr - min_year)
    rect = Rectangle((min_year, 0), width, 1, fc='silver', ec='white')
    ax2.add_patch(rect)
    ax2.set_xlim([min_year, max_year])
    ax2.set_ylim([0, 1])
    ax2.set_xticks(np.arange(min_year,max_year,5))
    ax2.set_yticks([])

    ax3 = fig.add_subplot(gs2[-3, 1:-2])
    width = (iter_number - min_iter) / (max_iter - min_iter)
    rect = Rectangle((date.year, 0), width, 1, fc='silver', ec='white')
    ax3.add_patch(rect)
    ax3.set_xlim([date.year, date.year + 1])
    ax3.set_ylim([0, 1])
    ax3.set_xticks(np.arange(date.year, date.year + 1, 1 / 12))
    ax3.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
    ax3.set_yticks([])

    output_path = os.path.join(output_dir, file_name)
    plt.savefig(output_path)
    plt.close(fig)

def compile_panels_to_movie(config_dir,config_name,field_name):
    pwd = os.getcwd()

    panels_dir = os.path.join(config_dir,'L1',config_name,'plots','results','panels', field_name)

    # get a list of the file names
    all_file_names = []
    for file_name in os.listdir(panels_dir):
        if file_name[0]!='.' and file_name[-4:]=='.png':
            all_file_names.append(file_name)
    all_file_names = sorted(all_file_names)

    # make a temporary dir where we will link all available images and go there
    panels_dir_tmp = os.path.join(config_dir, 'L1', config_name, 'plots', 'results', 'panels_tmp')
    os.mkdir(panels_dir_tmp)
    os.chdir(panels_dir_tmp)

    # link all of the images
    for ff in range(len(all_file_names)):
        os.system('ln -s ' + '../panels/'+field_name+'/'+all_file_names[ff]+' panel_'+'{:05d}'.format(ff)+'.png')

    output_name = config_name+'_'+field_name+'.mp4'

    os.system("ffmpeg -r 5 -i panel_%05d.png -vcodec mpeg4 -b 3M -y " + output_name)
    os.rename(output_name, os.path.join('..', output_name))

    os.chdir(pwd)
    os.system('rm -rf ' + panels_dir_tmp)

def create_variable_movies(config_dir, config_name, field_names, metadata_dict, remove_old, skip, print_level):

    # this can be implemented later
    # first need some sort of script to generate it as I did for the previous model
    add_background_imagery = False

    results_dir = os.path.join(config_dir, 'L1', config_name, 'results')

    for field_name in field_names:

        field_grid, all_iter_numbers = read_field_from_monthly_ncs(config_dir, config_name, results_dir, field_name)

        if 'results' not in os.listdir(os.path.join(config_dir,'L1',config_name,'plots',)):
            os.mkdir(os.path.join(config_dir,'L1',config_name,'plots','results'))

        if 'panels' not in os.listdir(os.path.join(config_dir,'L1',config_name,'plots','results')):
            os.mkdir(os.path.join(config_dir,'L1',config_name,'plots','results','panels'))

        output_dir = os.path.join(config_dir,'L1',config_name,'plots','results','panels')
        if field_name not in os.listdir(output_dir):
            os.mkdir(os.path.join(output_dir,field_name))
        output_dir = os.path.join(output_dir,field_name)

        all_dates = []
        for iter_number in all_iter_numbers:
            date = iter_number_to_date(iter_number[0])
            ymd_string = str(date.year)+'{:02d}'.format(date.month)+'{:02d}'.format(date.day)
            all_dates.append(ymd_string)

        ds = nc4.Dataset(os.path.join(config_dir, 'nc_grids', config_name + '_grid.nc'))
        depth = ds.variables['Depth'][:, :]
        ds.close()
        if field_name == 'Vorticity' or field_name == 'Vorticity_AW':
            depth = depth[1:, 1:]

        land_mask = (depth == 0).astype(int)
        land_mask = np.ma.masked_where(land_mask == 0, land_mask)

        aw_depth_mask = np.logical_and(depth > 0, depth < 257).astype(int)
        aw_depth_mask = np.ma.masked_where(aw_depth_mask == 0, aw_depth_mask)

        if remove_old:
            os.system('rm -rf '+output_dir+'/*')

        for i in range(len(all_iter_numbers)):

            output_file = config_name+'_'+field_name+'_'+all_dates[i]+'.png'

            if output_file not in os.listdir(output_dir):
                print('    - Creating plot for date '+str(all_dates[i]))

                iter_number = all_iter_numbers[i][0]

                create_panel_plot(output_dir, output_file, metadata_dict, field_name, field_grid, i, iter_number, land_mask, aw_depth_mask, depth, add_background_imagery)


        compile_panels_to_movie(config_dir, config_name, field_name)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-r", "--remove_old", action="store",
                        help="Choose whether to remove old files (1 is true, 0 is false).", dest="remove_old",
                        type=int, required=False, default = 0)

    parser.add_argument("-s", "--skip", action="store",
                        help="Choose how many panels to skip at a time.", dest="skip",
                        type=int, required=False, default=1)

    args = parser.parse_args()
    config_dir = args.config_dir
    remove_old = args.remove_old
    skip = args.skip

    if remove_old==0:
        remove_old = False
    else:
        remove_old = True

    create_variable_movies(config_dir,remove_old,skip)
   

