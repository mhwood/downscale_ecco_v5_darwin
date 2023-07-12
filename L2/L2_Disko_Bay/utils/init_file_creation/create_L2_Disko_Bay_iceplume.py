
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_L2_iceplume_files(config_dir, runoff_dir, termpicks_file, glacier_IDs, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    sys.path.insert(1, os.path.join(config_dir,'L2', 'utils','init_file_creation'))

    model_name = 'L2_Disko_Bay'

    print('    - Ensure these steps are completed to create the iceplume files')
    print('        1. Download Mankoff\'s runoff data from here: ')
    print('           https://dataverse.geus.dk/dataset.xhtml?persistentId=doi:10.22008/FK2/XKQVL7')
    print('        2. Make a shapefile of the domain using the utils/domain_shapefile.py script')
    print('        3. In GIS: manually identify near-front points as 1 (melting) or 3 (plume)')
    print('            - Save the new shapefile as ...')
    print('        4. Run the create_L2_iceplume_files script below')

    # import domain_shapefile as ds
    # ds.create_shp(config_dir, 'L2', model_name, print_level)

    print('DIDNT FINISH - READ CLOSELY TO FIX SUBGLACIAL DISCHARGE')
    years = np.arange(1994,2019).tolist()
    import create_L2_iceplume_files as cip
    cip.create_L2_iceplume_files(config_dir, model_name, runoff_dir, termpicks_file, glacier_IDs, years, print_level)








if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_L3_iceplume_files(config_dir)
   

