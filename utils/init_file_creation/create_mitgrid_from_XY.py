
import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import pyproj
import argparse


def create_mitgrid_matrices(domain_coords):
    mitgrid_matrices = dict()
    mitgrid_matrices['XG'] = domain_coords[2]
    mitgrid_matrices['YG'] = domain_coords[3]
    mitgrid_matrices['XC'] = domain_coords[0]
    mitgrid_matrices['YC'] = domain_coords[1]
    return(mitgrid_matrices)

def create_new_mitgrid(output_file,mitgrid_matrices,XG,YG,factor):
    mg_new, n_rows, n_cols = \
        sg.regrid.regrid(mitgrid_matrices=mitgrid_matrices, \
                         lon_subscale=factor, lat_subscale=factor, \
                         lon1=XG[0, 0], lat1=YG[0, 0],
                         lon2=XG[-1, -1], lat2=YG[-1, -1],
                         verbose=False, outfile = output_file)

    print('    - Output shape: ('+str(n_rows)+', '+str(n_cols)+')')

    sg.gridio.write_mitgridfile(output_file, mg_new, n_rows, n_cols)

def create_mitgrid_file(config_dir,model_level,model_name,XC,YC,XG,YG,print_level):

    domain_coords = [XC, YC, XG, YG]

    mitgrid_matrices = create_mitgrid_matrices(domain_coords)

    if print_level >=1:
        print('    - Generating mitgrid for the '+model_name+' domain')
        print('        - Simplegrid is slow, this may take minute or two')
    output_file = os.path.join(config_dir,model_level,model_name,'input',model_name+'.mitgrid')

    create_new_mitgrid(output_file, mitgrid_matrices, XG, YG, factor=1)
