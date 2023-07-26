
import os
import netCDF4 as nc4
import numpy as np
import argparse
import sys
import cmocean.cm as cm

def create_movies(config_dir, field_name, print_level):

    L1_model_name = 'L1_East_Pacific'

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils','plot_creation','results'))

    vel_fields = ['Speed','Vorticity']

    Chl_fields = ['Chl01','Chl02','Chl03','Chl04','Chl05']

    const_fields = ['DIC','FeT','NH4','NO2','NO3','PO4','SiO2']

    DO_fields = ['DOC','DOFe','DON','DOP']

    cx_fields = ['c01', 'c02', 'c03', 'c04', 'c05']

    PO_fields = ['POC','POFe','PON','POP','POSi']

    misc_fields = ['Alk','O2','PIC']

    TS_fields = ['Theta','Salt']

    etan_fields = ['EtaN']

    field_names = vel_fields + Chl_fields + const_fields + DO_fields + cx_fields + PO_fields + misc_fields + TS_fields + etan_fields

    field_names = [field_name]

    metadata_dict = {'EtaN': [0, 1, 'viridis', 'm', 'Surface Height Anomaly'],
                     'Theta': [10, 25, 'turbo', '$^{\circ}$C', 'Potential Temperature (Surface)'],  #
                     'Salt': [25, 35, cm.haline, 'psu', 'Practical Salinity (Surface)'],  #
                     'Uvel': [-1, 1, cm.balance, 'm/s', 'Zonal Velocity'],  #
                     'Vvel': [-1, 1, cm.balance, 'm/s', 'Meridional Velocity'],  #
                     'Speed': [0, 0.25, cm.tempo_r, 'm/s', 'Speed (Surface)'],
                     'Vorticity': [-0.15, 0.15, cm.balance, '$\zeta$/f', 'Vorticity (Surface)'],
                     'Vorticity_AW': [-0.4, 0.4, cm.balance, '$\zeta$/f', 'Subsurface Vorticity (257 m)'],
                     'DIC': [0, 1, cm.amp_r, '$\mu$M C', 'Dissolved Inorganic Carbon'],
                     'NO3': [0, 1, cm.dense, '$\mu$M N', 'Nitrate'],
                     'NO2': [0, 1, cm.dense, '$\mu$M N', 'Nitrite'],
                     'NH4': [0, 1, cm.dense, '$\mu$M N', 'Ammonium'],
                     'PO4': [0, 1, cm.matter, '$\mu$M P', 'Phosphate'],
                     'FeT': [0, 1, cm.turbid_r, '$\mu$M Fe', 'Terrestrial Iron'],
                     'SiO2': [0, 1, cm.gray, '$\mu$M Si', 'Silicon Dioxide'],
                     'DOC': [0, 50, 'turbo', '$\mu$M C', 'Dissolved Organic Carbon'],
                     'DON': [0, 1, cm.dense, '$\mu$M N', 'Dissolved Organic Nitrogen'],
                    'DOP': [0, 1, cm.matter, '$\mu$M P', 'Dissolved Organic Phosphorus'],
                    'DOFe': [0, 1, cm.turbid_r, '$\mu$M Fe', 'Dissolved Organic Iron'],
                    'POC': [0, 1, cm.amp_r, '$\mu$M C', 'Particulate Organic Carbon'],
                    'PON': [0, 1, cm.dense, '$\mu$M N', 'Particulate Organic Nitrogen'],
                    'POP': [0, 1, cm.matter, '$\mu$M P', 'Particulate Organic Phosphorus'],
                    'POFe': [0, 1, cm.turbid_r, '$\mu$M Fe', 'Particulate Organic Iron'],
                    'POSi': [0, 1, cm.gray, '$\mu$M Si', 'Particulate Organic Silicon'],
                    'PIC': [0, 1, cm.amp_r, '$\mu$M C', 'Particulate Inrganic Carbon'],
                    'Alk': [0, 1, cm.amp_r, 'meq/m$^3$', 'Alkalinity'],
                    'O2': [0, 1, cm.tempo_r, '$\mu$M O', 'Oxygen'],
                    'c01': [0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c02': [0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c03': [0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c04': [0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c05': [0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c06': [0, 1, cm.amp_r, '$\mu$M C', ''],
                    'c07': [0, 1, cm.amp_r, '$\mu$M C', ''],
                    'Chl01': [0, 1, 'turbo', 'mg/m$^3$', ''],
                    'Chl02': [0, 1, 'turbo', 'mg/m$^3$', ''],
                    'Chl03': [0, 1, 'turbo', 'mg/m$^3$', ''],
                    'Chl04': [0, 1, 'turbo', 'mg/m$^3$', ''],
                    'Chl05': [0, 1, 'turbo', 'mg/m$^3$', ''],
                    'Total_Chl': [0, 1, 'turbo', 'mg/m$^3$', '']}

    remove_old = False
    skip = False

    # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    import create_L1_variable_movies as pvm
    pvm.create_variable_movies(config_dir, L1_model_name, field_names, metadata_dict, remove_old, skip, print_level)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-f", "--field_name", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="field_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    field_name = args.field_name

    create_movies(config_dir, field_name, print_level=4)
   

