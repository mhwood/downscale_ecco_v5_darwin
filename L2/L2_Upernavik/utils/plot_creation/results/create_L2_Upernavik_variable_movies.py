
import os
import numpy as np
import netCDF4 as nc4
import argparse
import sys
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cmocean.cm as cm
import shutil
import argparse
from matplotlib.patches import Rectangle
import datetime as dt
from datetime import datetime, timedelta
from osgeo import gdal


def read_field_from_monthly_ncs(config_dir, year, results_dir, field_metadata, field_name, depth_index):

    if field_name == 'Total_Chl':
        var_set = ['Chl01','Chl02','Chl03','Chl04','Chl05']
    else:
        var_set = [field_name]

    grid_started = False

    year_months = []
    for file_name in os.listdir(os.path.join(results_dir, field_metadata[0], var_set[0])):
        if file_name[-3:]=='.nc' and file_name[0]!='.':

            year_month = file_name.split('_')[-1].split('.')[0]
            if year_month not in year_months and int(year_month[0:4])==year:# and int(year_month[4:6]) in [1]:
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

    grid = []
    all_iter_numbers = []
    for year_month in year_months:
        for file_name in os.listdir(os.path.join(results_dir,field_metadata[0], var_set[0])):
            # print(file_name)
            if file_name.split('_')[-1].split('.')[0]==year_month and file_name[-3:]=='.nc' and file_name[0]!='.':
                if field_name!='Total_Chl':
                    print('    - Reading from '+file_name)

                if field_name in ['Total_Chl']:
                    print('    - Reading from ' + file_name)
                    ds = nc4.Dataset(os.path.join(results_dir, field_metadata[0], var_set[0], file_name))
                    field = ds.variables['Chl01'][:, :, :, :]
                    ds.close()
                    field = field[:,depth_index,:,:]

                    print('    - Reading from ' + file_name.replace('Chl01','Chl02'))
                    ds = nc4.Dataset(os.path.join(results_dir, field_metadata[0],var_set[1], file_name.replace('Chl01','Chl02')))
                    chl2 = ds.variables['Chl02'][:, :, :, :]
                    ds.close()
                    field += chl2[:,depth_index,:,:]

                    print('    - Reading from ' + file_name.replace('Chl01', 'Chl03'))
                    ds = nc4.Dataset(os.path.join(results_dir, field_metadata[0],var_set[2], file_name.replace('Chl01', 'Chl03')))
                    chl2 = ds.variables['Chl03'][:, :, :, :]
                    ds.close()
                    field += chl2[:,depth_index,:,:]

                    print('    - Reading from ' + file_name.replace('Chl01', 'Chl04'))
                    ds = nc4.Dataset(os.path.join(results_dir, field_metadata[0],var_set[3], file_name.replace('Chl01', 'Chl04')))
                    chl2 = ds.variables['Chl04'][:, :, :, :]
                    ds.close()
                    field += chl2[:,depth_index,:,:]

                    print('    - Reading from ' + file_name.replace('Chl01', 'Chl05'))
                    ds = nc4.Dataset(os.path.join(results_dir, field_metadata[0],var_set[4], file_name.replace('Chl01', 'Chl05')))
                    chl2 = ds.variables['Chl05'][:, :, :, :]
                    iter_numbers = ds.variables['iterations'][:]
                    ds.close()
                    field += chl2[:,depth_index,:,:]

                else:
                    ds = nc4.Dataset(os.path.join(results_dir, field_metadata[0], field_name, file_name))
                    field = ds.variables[field_name][:, :, :]
                    iter_numbers = ds.variables['iterations'][:]
                    ds.close()
                    if len(np.shape(field))==4:
                        field = field[:,depth_index,:,:]
                    else:
                        field = field[:, :, :]

                if not grid_started:
                    grid_started = True
                    grid = field
                    all_iter_numbers = np.reshape(iter_numbers,(np.size(iter_numbers),1))
                else:
                    grid = np.concatenate([grid,field],axis=0)
                    all_iter_numbers = np.concatenate([all_iter_numbers,np.reshape(iter_numbers,(np.size(iter_numbers),1))],axis=0)

    return (grid, all_iter_numbers)

def read_background_imagery(file_path):

    ds = nc4.Dataset(file_path)
    R = ds.variables['band_1'][:,:]
    G = ds.variables['band_3'][:,:]
    B = ds.variables['band_4'][:,:]
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    ds.close()
    rows = np.shape(R)[0]
    cols = np.shape(R)[1]
    extents = [x[0], x[-1], y[0], y[-1]]


    R = R.reshape((rows,cols, 1))
    G = G.reshape((rows, cols, 1))
    B = B.reshape((rows, cols, 1))
    image = np.concatenate([B,G,R],axis=2)
    brightness_factor = 0.3 # 0 - 1
    image = (np.max(image)-np.max(image)*(brightness_factor))/(np.max(image)-np.min(image))*(image-np.min(image))+np.max(image)*(brightness_factor)
    image[image<0]=0
    image[image>1]=1

    # x_resolution = transform[1]
    # y_resolution = transform[5]
    return(image,extents)

def gef_field_metadata(field_name,depth_index):

    if field_name in ['Theta','Salt']:
        metadata_dict = {('Theta',39): ['daily_mean', 0.75, 1.5, 'turbo', '$^{\circ}$C', 'Potential Temperature'],
                         ('Salt',39): ['daily_mean', 34.1, 34.5, cm.haline, 'psu', 'Practical Salinity']}

        output = metadata_dict[(field_name, depth_index)]
    else:

        metadata_dict = {'EtaN': ['daily_snapshot', 0, 1, 'viridis', 'm', 'Surface Height Anomaly'],
                     'Theta': ['daily_mean', 10, 25, 'turbo', '$^{\circ}$C', 'Potential Temperature'],
                     'Theta_AW': ['daily_snapshot', -1, 2, 'turbo', '$^{\circ}$C', 'Potential Temperature'],#
                     'Theta_500_AW': ['daily_mean', -1, 2, 'turbo', '$^{\circ}$C', 'Potential Temperature (500 m)'],  #
                     'Salt': ['daily_mean', 25, 35, cm.haline, 'psu', 'Practical Salinity'],  #
                     'Salt_AW': ['daily_snapshot', 33.5, 34.5, cm.haline, 'psu', 'Practical Salinity'],  #
                     'Uvel': ['daily_mean', -1, 1, cm.balance, 'm/s', 'Zonal Velocity'],  #
                     'Vvel': ['daily_mean', -1, 1, cm.balance, 'm/s', 'Meridional Velocity'],  #
                     'Speed': ['daily_mean', 0, 0.25, cm.tempo_r, 'm/s', 'Speed'],
                     'Vorticity': ['daily_mean', -0.15, 0.15, cm.balance, '$\zeta$/f', 'Vorticity'],
                     'Vorticity_AW': ['daily_mean', -0.4, 0.4, cm.balance, '$\zeta$/f', 'Subsurface Vorticity'],
                     'DIC': ['daily_mean', 0, 1, cm.amp_r, '$\mu$M C', 'Dissolved Inorganic Carbon'],
                     'NO3': ['daily_mean', 3, 5, cm.dense, '$\mu$M N', 'Nitrate'],
                     'NO2': ['daily_mean', 0, 1, cm.dense, '$\mu$M N', 'Nitrite'],
                     'NH4': ['daily_mean', 0, 1, cm.dense, '$\mu$M N', 'Ammonium'],
                     'PO4': ['daily_mean', 0, 1, cm.matter, '$\mu$M P', 'Phosphate'],
                     'FeT': ['daily_mean', 0, 1, cm.turbid_r, '$\mu$M Fe', 'Terrestrial Iron'],
                     'SiO2': ['daily_mean', 0, 1, cm.gray, '$\mu$M Si', 'Silicon Dioxide'],
                     'DOC': ['daily_mean', 0, 50, 'turbo', '$\mu$M C', 'Dissolved Organic Carbon'],
                     'DON': ['daily_mean', 0, 1, cm.dense, '$\mu$M N', 'Dissolved Organic Nitrogen'],
                    'DOP': ['daily_mean', 0, 1, cm.matter, '$\mu$M P', 'Dissolved Organic Phosphorus'],
                    'DOFe': ['daily_mean', 0, 1, cm.turbid_r, '$\mu$M Fe', 'Dissolved Organic Iron'],
                    'POC': ['daily_mean', 0, 1, cm.amp_r, '$\mu$M C', 'Particulate Organic Carbon'],
                    'PON': ['daily_mean', 0, 1, cm.dense, '$\mu$M N', 'Particulate Organic Nitrogen'],
                    'POP': ['daily_mean', 0, 1, cm.matter, '$\mu$M P', 'Particulate Organic Phosphorus'],
                    'POFe': ['daily_mean', 0, 1, cm.turbid_r, '$\mu$M Fe', 'Particulate Organic Iron'],
                    'POSi': ['daily_mean', 0, 1, cm.gray, '$\mu$M Si', 'Particulate Organic Silicon'],
                    'PIC': ['daily_mean', 0, 1, cm.amp_r, '$\mu$M C', 'Particulate Inrganic Carbon'],
                    'Alk': ['daily_mean', 0, 1, cm.amp_r, 'meq/m$^3$', 'Alkalinity'],
                    'O2': ['daily_mean', 0, 1, cm.tempo_r, '$\mu$M O', 'Oxygen'],
                    'c01': ['daily_mean', 0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c02': ['daily_mean', 0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c03': ['daily_mean', 0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c04': ['daily_mean', 0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c05': ['daily_mean', 0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c06': ['daily_mean', 0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c07': ['daily_mean', 0, 1, cm.amp_r, '$\mu$M C', ''],
                    'Chl01': ['daily_mean', 0, 1, 'turbo', 'mg/m$^3$', 'Diatom Functional Group'],
                    'Chl02': ['daily_mean', 0, 1, 'turbo', 'mg/m$^3$', 'Large Eukaryote Functional Group'],
                    'Chl03': ['daily_mean', 0, 1, 'turbo', 'mg/m$^3$', 'Synechococcus Functional Group'],
                    'Chl04': ['daily_mean', 0, 1, 'turbo', 'mg/m$^3$', 'Low-light adapted Prochlorococcus Functional Group'],
                    'Chl05': ['daily_mean', 0, 1, 'turbo', 'mg/m$^3$', 'High-light adapted Prochlorococcus Functional Group'],
                     'Total_Chl': ['daily_mean', 0, 5, 'turbo', 'mg/m$^3$', 'Total Chlorophyll'],
                     'SIarea': ['daily_snapshot', 0, 1, cm.ice, '', '']}

        output = metadata_dict[field_name]

    return(output)

def iter_number_to_date(iter_number, seconds_per_iter):
    total_seconds = iter_number*seconds_per_iter
    date = datetime(1992,1,1) + timedelta(seconds=total_seconds)
    return(date)

def date_to_iter_number(date, seconds_per_iter):
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

def create_panel_plot(output_dir, year, experiment, file_name, config_name, field_metadata,
                      field_name, field_grid, i, iter_number, land_mask, aw_depth_mask, depth, depth_level, add_background_imagery):

    ############################################################################
    # get and generate some metadat information

    # metadata = metadata_dict[field_name]
    vmin = field_metadata[1]
    vmax = field_metadata[2]
    cmap = field_metadata[3]
    units = field_metadata[4]
    plot_title = field_metadata[5]

    if experiment in ['control','melange','baseline']:
        seconds_per_iter = 60
    else:
        seconds_per_iter = 30

    date = iter_number_to_date(iter_number, seconds_per_iter)
    year = date.year
    dec_yr = YMD_to_DecYr(date.year, date.month, date.day)
    min_iter = date_to_iter_number(datetime(date.year, 1, 1), seconds_per_iter)
    max_iter = date_to_iter_number(datetime(date.year + 1, 1, 1), seconds_per_iter)

    if add_background_imagery:
        file_path = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/configurations/' \
                    'downscale_darwin/L2/L2_Upernavik/plots_control/basemap/L2_Upernavik_MODIS_20220720_3413.nc'
        background_image, extents = read_background_imagery(file_path)

    ############################################################################
    # make the plot

    show_year_bar = False

    if show_year_bar:
        fig = plt.figure(figsize=(8.5, 7.5))
    else:
        fig = plt.figure(figsize=(7.5, 6.8))
    plt.style.use('dark_background')

    gs2 = GridSpec(17, 14, left=0.01, right=0.95, top = 0.95, bottom = 0.05, hspace=0.05)

    ax1 = fig.add_subplot(gs2[:-4, :-1])

    if add_background_imagery:
        rect = Rectangle((extents[0], extents[1]), extents[2], extents[3],
                         facecolor='white', edgecolor='white', zorder=1)
        ax1.add_patch(rect)
        ax1.imshow(background_image, extent=extents, alpha=0.7, zorder=1)

    plot_grid = np.copy(field_grid[i, :, :])
    if np.any(plot_grid!=0):
        print('        - range: ',np.min(plot_grid[plot_grid!=0]),np.max(plot_grid[plot_grid!=0]))
    if 'AW' in field_name:
        plot_grid = np.ma.masked_where(plot_grid==0,plot_grid)
        # land_mask = np.ma.masked_where(land_mask,plot_grid==0)
    else:
        plot_grid = np.ma.masked_where(depth <= 0, plot_grid)
    # plot_grid = np.ma.masked_where(plot_grid==0, plot_grid)
    C = ax1.imshow(plot_grid, extent=extents, origin='lower', vmin=vmin, vmax=vmax, cmap=cmap, zorder=1)

    # create a depth mask which is 1 where the plot grid is masked and where depth is greater than 0, and 0 elsewhere
    depth_mask = np.where(np.logical_or(depth>500, depth==0), 0, 1)
    depth_mask = np.ma.masked_where(depth_mask==0, depth_mask)
    # plot the depth mask in a black layer on top
    plt.imshow(depth_mask, extent=extents, origin='lower', vmin=0, vmax=1, cmap='Greys', zorder=2)

    # if 'AW' in field_name and add_background_imagery:
        # ax1.imshow(aw_depth_mask, origin='lower', vmin=0, vmax=1, cmap='Greys')
        # ax1.contour(depth,levels=[0.1],colors='w',linewidths=0.4)

    if not add_background_imagery:
        ax1.imshow(land_mask, origin='lower', vmin=0, vmax=2, cmap='Greys')

    ax1.set_xticks([])
    ax1.set_yticks([])

    cbar = plt.colorbar(C, fraction=0.026, pad=0.04)
    cbar.set_label(units, rotation=90)
    plt.title(plot_title +'(Depth: '+str(round(depth_level))+' m)')

    # ax1.text(5, 5, 'Timestep: ' + str(int(iter_number)), ha='left', va='bottom', color='k')

    min_year = 1992
    max_year = 2022

    if show_year_bar:
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
    n_days = 0
    xticks = []
    for month in range(1,13):
        if month in [1, 3, 5, 7, 8, 10, 12]:
            days_in_month = 31
        elif month in [4, 6, 9, 11]:
            days_in_month = 30
        else:
            if year % 4 == 0:
                days_in_month = 29
            else:
                days_in_month = 28
        xticks.append((n_days+days_in_month-15)/366 + date.year)
        n_days += days_in_month
        if month<12:
            plt.plot([n_days / 366 + date.year, n_days / 366 + date.year], [0, 1], color='white', linewidth=0.5)

    ax3.set_xticks(xticks)
    ax3.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
    ax3.set_yticks([])
    ax3.set_xlabel('Year '+str(date.year))

    output_path = os.path.join(output_dir, file_name)
    plt.savefig(output_path)
    plt.close(fig)

def compile_panels_to_movie(config_dir, plots_dir,config_name,field_name):
    pwd = os.getcwd()

    panels_dir = os.path.join(config_dir,'L2',config_name, plots_dir,'results','panels', field_name)

    # get a list of the file names
    all_file_names = []
    for file_name in os.listdir(panels_dir):
        if file_name[0]!='.' and file_name[-4:]=='.png':
            all_file_names.append(file_name)
    all_file_names = sorted(all_file_names)

    # make a temporary dir where we will link all available images and go there
    panels_dir_tmp = os.path.join(config_dir, 'L2', config_name,  plots_dir, 'results', 'panels_tmp')
    os.mkdir(panels_dir_tmp)
    os.chdir(panels_dir_tmp)

    # link all of the images
    for ff in range(len(all_file_names)):
        os.system('ln -s ' + '../panels/'+field_name+'/'+all_file_names[ff]+' panel_'+'{:05d}'.format(ff)+'.png')

    output_name = config_name+'_'+field_name+'.mp4'

    os.system("ffmpeg -r 10 -i panel_%05d.png -vcodec mpeg4 -b 3M -y " + output_name)
    os.rename(output_name, os.path.join('..', output_name))

    os.chdir(pwd)
    os.system('rm -rf ' + panels_dir_tmp)

def create_variable_movies(config_dir, plots_dir, year, experiment, field_name, depth_index,
                           field_metadata, remove_old, skip, print_level=4):

    # this can be implemented later
    # first need some sort of script to generate it as I did for the previous model
    add_background_imagery = True

    results_dir = os.path.join(config_dir, 'L2', 'L2_Upernavik', 'results_'+experiment)

    # if 'AW' in field_name:
    #     depth_index = 31
    # elif 'Chl' in field_name:
    #     depth_index = 10
    # else:
    #     depth_index = 0

    if experiment in ['control','melange','baseline']:
        seconds_per_iter = 60
    else:
        seconds_per_iter = 30

    field_grid, all_iter_numbers = read_field_from_monthly_ncs(config_dir, year, results_dir,
                                                               field_metadata, field_name, depth_index)
    print('    field grid shape:', np.shape(field_grid))
    print('    field metadata:', field_metadata)

    if 'results' not in os.listdir(os.path.join(config_dir,'L2','L2_Upernavik',plots_dir)):
        os.mkdir(os.path.join(config_dir,'L2','L2_Upernavik',plots_dir,'results'))

    if 'panels' not in os.listdir(os.path.join(config_dir,'L2','L2_Upernavik',plots_dir,'results')):
        os.mkdir(os.path.join(config_dir,'L2','L2_Upernavik',plots_dir,'results','panels'))

    output_dir = os.path.join(config_dir,'L2','L2_Upernavik',plots_dir,'results','panels')
    if field_name not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,field_name))
    output_dir = os.path.join(output_dir,field_name)

    all_dates = []
    for iter_number in all_iter_numbers:

        date = iter_number_to_date(iter_number[0], seconds_per_iter)
        ymd_string = str(date.year)+'{:02d}'.format(date.month)+'{:02d}'.format(date.day)
        all_dates.append(ymd_string)

    ds = nc4.Dataset(os.path.join(config_dir, 'nc_grids', 'L2_Upernavik' + '_grid.nc'))
    depth = ds.variables['Depth'][:, :]
    drF = ds.variables['drF'][:]
    Z = np.cumsum(drF)
    ds.close()
    if field_name == 'Vorticity' or field_name == 'Vorticity_AW':
        depth = depth[1:, 1:]
        depth = depth[1:, 1:]

    land_mask = (depth == 0).astype(int)
    land_mask = np.ma.masked_where(land_mask == 0, land_mask)

    aw_depth_mask = np.logical_and(depth > 0, depth < 257).astype(int)
    aw_depth_mask = np.ma.masked_where(aw_depth_mask == 0, aw_depth_mask)

    depth_level = Z[depth_index]

    if remove_old:
        os.system('rm -rf '+output_dir+'/*')

    for i in range(len(all_iter_numbers)):

        output_file = 'L2_Upernavik'+'_'+field_name+'_'+all_dates[i]+'.png'

        if output_file not in os.listdir(output_dir):
            print('    - Creating plot for date '+str(all_dates[i]))

            iter_number = all_iter_numbers[i][0]

            create_panel_plot(output_dir, year, experiment, output_file, 'L2_Upernavik', field_metadata,
                              field_name, field_grid, i, iter_number, land_mask, aw_depth_mask, depth, depth_level, add_background_imagery)

    compile_panels_to_movie(config_dir,  plots_dir, 'L2_Upernavik', field_name)


config_dir = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'darwin3/configurations/downscale_darwin'

experiment = 'baseline'
year = 2021

plots_dir = os.path.join(config_dir,'L2', 'L2_Upernavik', 'plots_'+experiment)

remove_old = False
skip = False

field_name = 'Salt'
depth_index = 39

field_metadata = gef_field_metadata(field_name,depth_index)

create_variable_movies(config_dir, plots_dir, year, experiment, field_name, depth_index,
                       field_metadata, remove_old, skip, print_level=4)

