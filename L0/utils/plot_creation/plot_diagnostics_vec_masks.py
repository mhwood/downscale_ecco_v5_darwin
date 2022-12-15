import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import argparse

def read_bathy_to_tiles(bathy_file,llc):

    tiles = []

    dimFacets = [llc, llc*3, llc, llc*3, llc, llc, llc*3, llc, llc*3, llc]

    bathy = np.fromfile(bathy_file,dtype='>f4')
    # bathy[bathy>0]=0
    # bathy[bathy<0]=-1

    total_points_read = 0
    for i in range(int(len(dimFacets)/2)):
        dim_1 = dimFacets[2*i]
        dim_2 = dimFacets[2*i+1]
        face = bathy[total_points_read:total_points_read+dim_1*dim_2]

        face = np.reshape(face,(dim_2,dim_1))

        if dim_1>llc:
            n_tiles = int(dim_1 / llc)
            for j in range(n_tiles):
                tile = face[:,llc * j:llc * (j + 1)]
                tiles.append(tile)
        elif dim_2>llc:
            n_tiles = int(dim_2 / llc)
            for j in range(n_tiles):
                tile = face[llc * j:llc * (j + 1),:]
                tiles.append(tile)
        else:
            tile = face
            tiles.append(tile)

        total_points_read += np.size(face)

    return(tiles)

def read_proc_counter_to_tiles(tile_proc_file,tiles):
    f = open(tile_proc_file)
    lines = f.read()
    f.close()
    lines=lines.split('\n')
    lines.pop(0)

    # total_proc_number,proc_number,face,tile,subtile_row,subtile_col
    proc_dict = {}
    for tile_number in range(len(tiles)):
        tile_set = []
        for line in lines:
            line = line.split(',')
            if int(line[3])==tile_number+1:
                tile_set.append([int(line[1]),int(line[4]),int(line[5])])
        proc_dict[tile_number+1] = tile_set
    return(proc_dict)

def read_masks_to_tiles(mask_dir,llc,mask_names):
    mask_sets = []
    for mask_name in mask_names:
        mask_file = os.path.join(mask_dir,mask_name)
        mask_tiles = read_bathy_to_tiles(mask_file,llc)
        mask_sets.append(mask_tiles)
    return(mask_sets)


def create_plot(output_file,tiles,llc,mask_names,mask_set):
    tile_to_subplot = {1: [5, 1], 2: [4, 1], 3: [3, 1],
                       4: [5, 2], 5: [4, 2], 6: [3, 2],
                       7: [2, 2],
                       8: [2, 3], 9: [2, 4], 10: [2, 5],
                       11: [1, 3], 12: [1, 4], 13: [1, 5]}

    fig = plt.figure(figsize=(15, 15))
    plt.style.use('dark_background')

    vmin=0
    vmax=0
    for tile in tiles:
        if np.min(tile)<vmin:
            vmin=np.min(tile)
    vmin=0.8*vmin

    vmin = -500
    vmax = 0

    mask_colors = ['red','orange','green','purple']

    proc_counter = 0

    for tile_number in range(len(tiles)):
        row = tile_to_subplot[tile_number+1][0]
        col = tile_to_subplot[tile_number+1][1]
        counter = (row - 1) * 5 + col
        plt.subplot(5, 5, counter)

        bathy_subset = tiles[tile_number]
        plt.imshow(bathy_subset, origin='lower', cmap='Blues_r',vmin=vmin,vmax=vmax)

        for mi in range(len(mask_names)):
            mask_tiles = mask_set[mi]
            mask_subset = mask_tiles[tile_number]
            if np.any(mask_subset>0):
                plot_row,plot_col = np.where(mask_subset>0)
                for ri in range(len(plot_row)):
                    plt.plot(plot_col[ri],plot_row[ri],'.',color=mask_colors[mi],markersize=2)

        for proc_row in range(3):
            for proc_col in range(3):
                proc_counter+=1
                plt.text(proc_col*90+45,proc_row*90+45,str(proc_counter),fontsize=14,color='k',va='center',ha='center')
                plt.text(proc_col * 90 + 45, proc_row * 90 + 45, str(proc_counter),fontsize=11,va='center',ha='center')

        plt.plot([90,90],[0,270],'k-')
        plt.plot([180, 180], [0, 270],'k-')
        plt.plot([0, 270],[90, 90],  'k-')
        plt.plot([0, 270],[180, 180],  'k-')

        plt.gca().set_xlim([0, llc])
        plt.gca().set_ylim([0, llc])

    plt.subplot(5, 5, 19)
    for mi in range(len(mask_names)):
        plt.plot([0,1],[mi,mi],'-',color=mask_colors[mi])
        plt.text(1.1,mi,mask_names[mi],fontsize=12,ha='left',va='center')

    plt.gca().set_xlim([-0.1,3])
    plt.gca().set_ylim([-1,len(mask_names)])
    plt.axis('off')

    plt.savefig(output_file,bbox_inches='tight')
    plt.close(fig)

########################################################################################################################

def plot_dv_masks(ecco_path,model_names):
    llc = 270

    if 'plots' not in os.listdir(os.path.join('..')):
        os.mkdir(os.path.join('..','plots'))

    bathy_file = os.path.join(ecco_path,'LLC270_Files','input_init','bathy_llc270')

    tiles = read_bathy_to_tiles(bathy_file,llc)
    for tile in tiles:
        tile[tile>0]=0

    mask_dir = os.path.join('..', 'input','dv')

    for model_name in model_names:

        mask_names = []
        for boundary in ['north','south','east','west']:
            file_name = model_name+'_'+boundary+'.bin'
            if file_name in os.listdir(mask_dir):
                mask_names.append(model_name+'_'+boundary+'.bin')

        mask_set = read_masks_to_tiles(mask_dir,llc,mask_names)
        output_file = os.path.join('..', 'plots', model_name+'_boundary_dv_masks.png')
        create_plot(output_file,tiles,llc,mask_names,mask_set)

        mask_set = read_masks_to_tiles(mask_dir, llc, [model_name+'_surf.bin'])
        output_file = os.path.join('..', 'plots', model_name + '_surface_dv_mask.png')
        create_plot(output_file, tiles, llc, [model_name+'_surface_mask'], mask_set)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=True)

    parser.add_argument("-m", "--model_names", action="store", help="List of models to plot.",
                        dest="model_names", type=str, nargs='+', required=True)

    args = parser.parse_args()
    ecco_path = args.ecco_path
    model_names = args.model_names

    plot_dv_masks(ecco_path,model_names)
