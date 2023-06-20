
import os
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import cmocean.cm as cm
import datetime


def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_boundary_condition(config_dir, L1_model_name, boundaries, field_name, level='L1', balanced = True):

    dec_yr = []

    # make the dec yr array
    for year in range(1992,2023):
        for month in range(1,13):
            if month in [1, 3, 5, 7, 8, 10, 12]:
                nDays = 31
            elif month in [4, 6, 9, 11]:
                nDays = 30
            else:
                if year % 4 == 0:
                    nDays = 29
                else:
                    nDays = 28
            for d in range(1,nDays+1):
                dec_yr.append(YMD_to_DecYr(year,month,d))

    time = np.array(dec_yr)

    if level=='L1' and balanced:
        file_path = os.path.join(config_dir,'L1',L1_model_name,'input','obcs',
                                     'timeseries','L1_timeseries','balanced',L1_model_name+'_'+field_name+'_boundary_flux.nc')
    if level=='L1' and not balanced:
        file_path = os.path.join(config_dir,'L1',L1_model_name,'input','obcs',
                                     'timeseries','L1_timeseries','unbalanced',L1_model_name+'_'+field_name+'_boundary_flux.nc')
    if level=='L0':
        file_path = os.path.join(config_dir,'L1',L1_model_name,'input','obcs',
                                     'timeseries','L0_timeseries','L0_'+L1_model_name[3:]+'_'+field_name+'_boundary_flux.nc')

    all_timeseries = []
    for boundary in boundaries:
        if os.path.exists(file_path):
            ds = nc4.Dataset(file_path)
            grp = ds.groups[boundary]
            timeseries = grp.variables['integrated_timeseries'][:]
            ds.close()
            time = time[:len(timeseries)]
        else:
            timeseries = []
        all_timeseries.append(timeseries)

    return(time,all_timeseries)

def calculate_integrated_flux_anomaly(time,timeseries):

    anomaly_timeseries = timeseries - np.mean(timeseries)

    timestep = 24*60*60 # 1 day

    integrated_anomaly = np.zeros_like(timeseries)
    for i in range(len(timeseries)-1):
        integrated_anomaly[i+1] = integrated_anomaly[i] + anomaly_timeseries[i]*timestep

    return(integrated_anomaly)

def tracer_number_to_metadata(tracer_number):
    tracer_dict = {1:['DIC','$\mu$M C',cm.amp_r],
                   2:['NO$_3$','$\mu$M N',cm.dense],
                   3:['NO$_2$','$\mu$M N',cm.dense],
                   4:['NH$_4$','$\mu$M N',cm.dense],
                   5:['PO$_4$','$\mu$M P',cm.matter],
                   6:['Fe$_T$','$\mu$M Fe',cm.turbid_r],
                   7:['SiO$_2$','$\mu$M Si',cm.gray],
                   8:['DOC','$\mu$M C',cm.amp_r],
                   9:['DON','$\mu$M N',cm.dense],
                   10:['DOP','$\mu$M P',cm.matter],
                   11:['DOFe','$\mu$M Fe',cm.turbid_r],
                   12:['POC','$\mu$M C',cm.amp_r],
                   13:['PON','$\mu$M N',cm.dense],
                   14:['POP','$\mu$M P',cm.matter],
                   15:['POFe','$\mu$M Fe',cm.turbid_r],
                   16:['POSi','$\mu$M Si',cm.gray],
                   17:['PIC','$\mu$M C',cm.amp_r],
                   18:['Alk','$\mu$M meq/m$^3$',cm.amp_r],
                   19:['O$_2$','$\mu$M O',cm.tempo_r],
                   20:['c01','$\mu$M C',cm.amp_r],
                   21:['c02','$\mu$M C',cm.amp_r],
                   22:['c03','$\mu$M C',cm.amp_r],
                   23:['c04','$\mu$M C',cm.amp_r],
                   24:['c05','$\mu$M C',cm.amp_r],
                   25:['c06','$\mu$M C',cm.amp_r],
                   26:['c07','$\mu$M C',cm.amp_r],
                   27:['Chl01','mg/m$^3$',cm.algae],
                   28:['Chl02','mg/m$^3$',cm.algae],
                   29:['Chl03','mg/m$^3$',cm.algae],
                   30:['Chl04','mg/m$^3$',cm.algae],
                   31:['Chl05','mg/m$^3$',cm.algae]}
    return(tracer_dict[tracer_number])

def plot_BC_timeseries_comparison(config_dir, L1_model_name, boundaries, velocity_only = False):

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'plot_creation','init_files'))

    if 'plots' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name)):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'plots'))
    if 'init_files' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name,'plots')):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'plots','init_files'))
    if 'BCs_L0_vs_L1_timeseries' not in os.listdir(os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files')):
        os.mkdir(os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files','BCs_L0_vs_L1_timeseries'))

    if velocity_only:
        var_names = ['VEL']
    else:
        if L1_model_name in ['L1_W_Greenland', 'L1_mac_delta']:
            var_names = ['THETA', 'SALT']#, 'VEL', 'UICE', 'VICE', 'HSNOW', 'HEFF', 'AREA']
        else:
            var_names = ['THETA', 'SALT', 'VEL']
        for p in range(1,32):
            var_names.append('PTRACE'+'{:02d}'.format(p))
        # var_names = ['PTRACE07']

    net_started = False

    for var_name in var_names:
        print('    - Plotting the '+var_name +' boundary conditions')

        subplot_count = 1

        fig = plt.figure(figsize=(14,9))
        plt.style.use('dark_background')

        time, all_L1_unbalanced_timeseries = read_boundary_condition(config_dir, L1_model_name, boundaries, var_name,balanced=False)
        time, all_L1_balanced_timeseries = read_boundary_condition(config_dir, L1_model_name, boundaries,var_name, balanced=True)
        time, all_L0_timeseries = read_boundary_condition(config_dir, L1_model_name, boundaries, var_name, level='L0')

        all_L1_unbalanced_anomalies = []
        all_L1_balanced_anomalies = []
        all_L0_anomalies = []
        for i in range(len(all_L1_unbalanced_timeseries)):
            L1_unbalanced_flux_anomaly = calculate_integrated_flux_anomaly(time,all_L1_unbalanced_timeseries[i])
            L1_balanced_flux_anomaly = calculate_integrated_flux_anomaly(time, all_L1_balanced_timeseries[i])
            L0_flux_anomaly = calculate_integrated_flux_anomaly(time, all_L0_timeseries[i])

            all_L1_unbalanced_anomalies.append(L1_unbalanced_flux_anomaly)
            all_L1_balanced_anomalies.append(L1_balanced_flux_anomaly)
            all_L0_anomalies.append(L0_flux_anomaly)

        #########################################################################################################

        if 'north' in boundaries:
            print('    - Plotting the timeseries on the north boundary')
            plt.subplot(len(boundaries), 2, subplot_count)
            north_L0_timeseries = all_L0_timeseries[boundaries.index('north')]
            north_L1_unbalanced_timeseries = all_L1_unbalanced_timeseries[boundaries.index('north')]
            north_L1_balanced_timeseries = all_L1_balanced_timeseries[boundaries.index('north')]
            plt.plot(time,north_L0_timeseries,label='L0', alpha=0.5)
            plt.plot(time, north_L1_unbalanced_timeseries, label='L1 unbalanced', alpha=0.5)
            if len(north_L1_balanced_timeseries)>0:
                plt.plot(time, north_L1_balanced_timeseries, label='L1 balanced', alpha=0.5)
            plt.grid(linestyle='--',linewidth=0.5,alpha=0.4)
            plt.ylabel('North')
            plt.title('Cumulative Flux: ' + str(var_name))
            if 'east' in boundaries or 'west' in boundaries or 'south' in boundaries:
                plt.gca().set_xticklabels([])
            plt.legend()
            subplot_count +=1

            plt.subplot(len(boundaries), 2, subplot_count)
            north_L0_anomalies = all_L0_anomalies[boundaries.index('north')]
            north_L1_unbalanced_anomalies = all_L1_unbalanced_anomalies[boundaries.index('north')]
            north_L1_balanced_anomalies = all_L1_balanced_anomalies[boundaries.index('north')]
            plt.plot(time, north_L0_anomalies, label='L0', alpha=0.5)
            plt.plot(time, north_L1_unbalanced_anomalies, label='L1 unbalanced', alpha=0.5)
            if len(north_L1_balanced_anomalies) > 0:
                plt.plot(time, north_L1_balanced_anomalies, label='L1 balanced', alpha=0.5)
            plt.grid(linestyle='--', linewidth=0.5, alpha=0.4)
            plt.ylabel('North')
            plt.title('Integrated Flux Anomaly: '+str(var_name))
            if 'east' in boundaries or 'west' in boundaries or 'south' in boundaries:
                plt.gca().set_xticklabels([])
            plt.legend()
            subplot_count += 1

        if 'south' in boundaries:
            print('    - Plotting the timeseries on the south boundary')
            plt.subplot(len(boundaries),2,subplot_count)
            south_L0_timeseries = all_L0_timeseries[boundaries.index('south')]
            south_L1_unbalanced_timeseries = all_L1_unbalanced_timeseries[boundaries.index('south')]
            south_L1_balanced_timeseries = all_L1_balanced_timeseries[boundaries.index('south')]
            plt.plot(time, south_L0_timeseries, label='L0', alpha=0.5)
            plt.plot(time, south_L1_unbalanced_timeseries, label='L1 unbalanced', alpha=0.5)
            if len(south_L1_balanced_timeseries)>0:
                plt.plot(time, south_L1_balanced_timeseries, label='L1 balanced', alpha=0.5)
            plt.grid(linestyle='--', linewidth=0.5, alpha=0.4)
            plt.ylabel('South')
            if 'east' in boundaries or 'west' in boundaries:
                plt.gca().set_xticklabels([])
            subplot_count +=1

            plt.subplot(len(boundaries), 2, subplot_count)
            south_L0_anomalies = all_L0_anomalies[boundaries.index('south')]
            south_L1_unbalanced_anomalies = all_L1_unbalanced_anomalies[boundaries.index('south')]
            south_L1_balanced_anomalies = all_L1_balanced_anomalies[boundaries.index('south')]
            plt.plot(time, south_L0_anomalies, label='L0', alpha=0.5)
            plt.plot(time, south_L1_unbalanced_anomalies, label='L1 unbalanced', alpha=0.5)
            if len(south_L1_balanced_anomalies) > 0:
                plt.plot(time, south_L1_balanced_anomalies, label='L1 balanced', alpha=0.5)
            plt.grid(linestyle='--', linewidth=0.5, alpha=0.4)
            plt.ylabel('South')
            if 'east' in boundaries or 'west' in boundaries:
                plt.gca().set_xticklabels([])
            subplot_count += 1

        if 'west' in boundaries:
            print('    - Plotting the timeseries on the west boundary')
            plt.subplot(len(boundaries),2,subplot_count)
            west_L0_timeseries = all_L0_timeseries[boundaries.index('west')]
            west_L1_unbalanced_timeseries = all_L1_unbalanced_timeseries[boundaries.index('west')]
            west_L1_balanced_timeseries = all_L1_balanced_timeseries[boundaries.index('west')]
            plt.plot(time, west_L0_timeseries, label='L0', alpha=0.5)
            plt.plot(time, west_L1_unbalanced_timeseries, label='L1 unbalanced', alpha=0.5)
            if len(west_L1_balanced_timeseries)>0:
                plt.plot(time, west_L1_balanced_timeseries, label='L1 balanced', alpha=0.5)
            plt.grid(linestyle='--', linewidth=0.5, alpha=0.4)
            plt.ylabel('West')
            if 'east' in boundaries:
                plt.gca().set_xticklabels([])
            subplot_count +=1

            plt.subplot(len(boundaries), 2, subplot_count)
            west_L0_anomalies = all_L0_anomalies[boundaries.index('west')]
            west_L1_unbalanced_anomalies = all_L1_unbalanced_anomalies[boundaries.index('west')]
            west_L1_balanced_anomalies = all_L1_balanced_anomalies[boundaries.index('west')]
            plt.plot(time, west_L0_anomalies, label='L0', alpha=0.5)
            plt.plot(time, west_L1_unbalanced_anomalies, label='L1 unbalanced', alpha=0.5)
            if len(west_L1_balanced_anomalies) > 0:
                plt.plot(time, west_L1_balanced_anomalies, label='L1 balanced', alpha=0.5)
            plt.grid(linestyle='--', linewidth=0.5, alpha=0.4)
            plt.ylabel('West')
            if 'east' in boundaries:
                plt.gca().set_xticklabels([])
            subplot_count += 1

        if 'east' in boundaries:
            print('    - Plotting the timeseries on the east boundary')
            plt.subplot(len(boundaries),2,subplot_count)
            east_L0_timeseries = all_L0_timeseries[boundaries.index('east')]
            east_L1_unbalanced_timeseries = all_L1_unbalanced_timeseries[boundaries.index('east')]
            east_L1_balanced_timeseries = all_L1_balanced_timeseries[boundaries.index('east')]
            plt.plot(time, east_L0_timeseries, label='L0', alpha=0.5)
            plt.plot(time, east_L1_unbalanced_timeseries, label='L1 unbalanced', alpha=0.5)
            if len(east_L1_balanced_timeseries)>0:
                plt.plot(time, east_L1_balanced_timeseries, label='L1 balanced', alpha=0.5)
            plt.grid(linestyle='--', linewidth=0.5, alpha=0.4)
            plt.ylabel('East')
            subplot_count += 1

            plt.subplot(len(boundaries), 2, subplot_count)
            east_L0_anomalies = all_L0_anomalies[boundaries.index('east')]
            east_L1_unbalanced_anomalies = all_L1_unbalanced_anomalies[boundaries.index('east')]
            east_L1_balanced_anomalies = all_L1_balanced_anomalies[boundaries.index('east')]
            plt.plot(time, east_L0_anomalies, label='L0', alpha=0.5)
            plt.plot(time, east_L1_unbalanced_anomalies, label='L1 unbalanced', alpha=0.5)
            if len(east_L1_balanced_anomalies) > 0:
                plt.plot(time, east_L1_balanced_anomalies, label='L1 balanced', alpha=0.5)
            plt.grid(linestyle='--', linewidth=0.5, alpha=0.4)
            plt.ylabel('East')

        plt.savefig(os.path.join(config_dir,'L1',L1_model_name,'plots', 'init_files','BCs_L0_vs_L1_timeseries',L1_model_name+'_BCs_L0_vs_L1_timeseries_'+var_name+'.png'))
        plt.close(fig)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_BC_timeseries_comparison(config_dir)
   

