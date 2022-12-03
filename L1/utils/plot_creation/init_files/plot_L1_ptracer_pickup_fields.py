
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from MITgcmutils import mds
import cmocean.cm as cm
import argparse
import sys



def read_ptracer_pickup_file(input_dir,pickup_iteration):
    print(' Reading from ' + os.path.join(input_dir, 'pickup_ptracers.' + '{:010d}'.format(pickup_iteration)))
    global_data, _, global_metadata = mds.rdmds(os.path.join(input_dir, 'pickup_ptracers.'+'{:010d}'.format(pickup_iteration)), returnmeta=True)

    var_names = []
    var_grids = []

    for vn in range(len(global_metadata['fldlist'])):
        var_grid = global_data[vn,0,:,:]
        var_grids.append(var_grid)
        var_names.append(global_metadata['fldlist'][vn].strip())

    return(var_names,var_grids,global_metadata)

def tracer_number_to_metadata(tracer_number):
    tracer_dict = {1:['DIC',cm.amp_r],
                   2:['NO$_3$',cm.dense],
                   3:['NO$_2$',cm.dense],
                   4:['NH$_4$',cm.dense],
                   5:['PO$_4$',cm.matter],
                   6:['Fe$_T$',cm.turbid_r],
                   7:['SiO$_2$',cm.gray],
                   8:['DOC',cm.amp_r],
                   9:['DON',cm.dense],
                   10:['DOP',cm.matter],
                   11:['DOFe',cm.turbid_r],
                   12:['POC',cm.amp_r],
                   13:['PON',cm.dense],
                   14:['POP',cm.matter],
                   15:['POFe',cm.turbid_r],
                   16:['POSi',cm.gray],
                   17:['PIC',cm.amp_r],
                   18:['Alk',cm.amp_r],
                   19:['O$_2$',cm.tempo_r],
                   20:['c01',cm.amp_r],
                   21:['c02',cm.amp_r],
                   22:['c03',cm.amp_r],
                   23:['c04',cm.amp_r],
                   24:['c05',cm.amp_r],
                   25:['c06',cm.amp_r],
                   26:['c07',cm.amp_r],
                   27:['Chl01',cm.algae],
                   28:['Chl02',cm.algae],
                   29:['Chl03',cm.algae],
                   30:['Chl04',cm.algae],
                   31:['Chl05',cm.algae]}
    return(tracer_dict[tracer_number])

def create_ptracer_pickup_plot(config_dir, L1_model_name, pickup_iteration):

    input_dir = os.path.join(config_dir,'L1', L1_model_name, 'input')

    var_names, var_grids, global_metadata = read_ptracer_pickup_file(input_dir,pickup_iteration)

    fig = plt.figure(figsize=(24, 12))
    plt.style.use("dark_background")

    counter = 1

    for ff in range(len(var_names)):

        plt.subplot(4, 8, counter)

        print(' Plotting ' + var_names[ff])
        var_grid_subset = var_grids[ff]

        metadata = tracer_number_to_metadata(ff+1)
        tracer_name = metadata[0]
        cmap = metadata[1]

        vmin = np.min(var_grid_subset[var_grid_subset != 0])
        vmax = np.max(var_grid_subset[var_grid_subset != 0])
        if vmin == vmax:
            vmin = -0.1
            vmax = 0.1

        C = plt.imshow(var_grid_subset, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        plt.colorbar(C, fraction=0.04, pad=0.04)
        plt.title(tracer_name)

        counter += 1

    output_file = os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files',
                               L1_model_name + '_ptracer_pickup_fields.png')
    plt.savefig(output_file, bbox_inches='tight')
    plt.close(fig)


########################################################################################################################
# User inputs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_ptracer_pickup_plot(config_dir)

