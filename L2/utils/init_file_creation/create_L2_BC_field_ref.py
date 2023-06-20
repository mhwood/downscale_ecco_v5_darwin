
import os
from datetime import datetime
import numpy as np
import netCDF4 as nc4
import argparse
import sys
import ast

def find_dv_files_to_read(config_dir, parent_model_level, parent_model_name, boundary, var_name, averaging_period, seconds_per_iter, points_per_output, Nr):

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import time_functions as tf

    iters_per_output = averaging_period/seconds_per_iter
    print(iters_per_output,averaging_period,seconds_per_iter)

    # loop through the dv files and mask a list of how many fields each file has
    file_names = []
    file_iters = []
    n_files_in_files = []
    dv_dir = os.path.join(config_dir, parent_model_level, parent_model_name, 'run', 'dv','L3_'+boundary)
    for file_name in os.listdir(dv_dir):
        if 'L3_' in file_name and var_name in file_name and boundary in file_name and file_name[0]!='.':
            file_iter = int(file_name.split('.')[-2])
            iters_per_file = int(np.size(np.fromfile(os.path.join(dv_dir,file_name),'>f4'))/(points_per_output*Nr))
            # iter_midpoints = tf.dv_file_name_iter_to_iter_midpoints(file_iter, iters_per_output, iters_per_file)
            file_iters.append(file_iter)
            file_names.append(file_name)
            n_files_in_files.append(iters_per_file)

    # sort the lists in case they got out of whack
    sorted_indices = sorted(range(len(file_iters)), key=lambda k: file_iters[k])
    file_names_sorted = []
    file_iters_sorted = []
    n_files_in_files_sorted = []
    for index in sorted_indices:
        file_names_sorted.append(file_names[index])
        file_iters_sorted.append(file_iters[index])
        n_files_in_files_sorted.append(n_files_in_files[index])

    # make dicts of the index and corresponding iters (w.r.t the parent model) which cen be read from these files
    iter_midpoint_dict = {}
    for i in range(len(file_iters_sorted)):
        file_iter = file_iters_sorted[i]
        file_name = file_names_sorted[i]
        iters_per_file = n_files_in_files_sorted[i]

        file_endpoints = np.arange(file_iter - 1, file_iter - 1 + (iters_per_file + 1) * iters_per_output, iters_per_output)
        file_midpoints = file_endpoints[:-1] + np.diff(file_endpoints) / 2

        # print(file_name,file_iter,iters_per_file)
        # iter_midpoints = tf.dv_file_name_iter_to_iter_midpoints(file_iter, iters_per_output, iters_per_file)
        iter_midpoint_dict[file_name] = file_midpoints

    return(file_names, iter_midpoint_dict)

def create_destination_file_list(config_dir, var_name, file_names, iter_midpoint_dict, averaging_period, seconds_per_iter, print_level):
    sys.path.insert(1, os.path.join(config_dir, 'utils'))
    import time_functions as tf

    start_seconds = 0
    iters_per_output = averaging_period / seconds_per_iter

    # create a list of daily bounds
    date_bounds = []
    for year in range(1992,2022):
        for month in range(1, 13):
            if month in [1, 3, 5, 7, 8, 10, 12]:
                nDays = 31
            elif month in [4, 6, 9, 11]:
                nDays = 30
            else:
                if year % 4 == 0:
                    nDays = 29
                else:
                    nDays = 28
            start_day = 1
            if year==1992:
                start_day = 2
            for day in range(start_day, nDays):
                date_bounds.append([datetime(year, month, day), datetime(year, month, day + 1)])
            if month != 12:
                date_bounds.append([datetime(year, month, day + 1), datetime(year, month + 1, 1)])
        date_bounds.append([datetime(year, 12, 31), datetime(year + 1, 1, 1)])

    dest_files = []
    dest_file_iter_bounds = {}
    # convert these to iteration numbers
    for date_set in date_bounds:
        date_0 = (date_set[0] - datetime(1992, 1, 1)).total_seconds()
        date_1 = (date_set[1] - datetime(1992, 1, 1)).total_seconds()
        iter_0 = tf.elapsed_seconds_to_iters(start_seconds, seconds_per_iter, date_0)
        iter_1 = tf.elapsed_seconds_to_iters(start_seconds, seconds_per_iter, date_1)
        dest_file = 'L3_BC_'+str(var_name)+'.'+str(date_set[0].year)+'{:02d}'.format(date_set[0].month)+'{:02d}'.format(date_set[0].day)+'.bin'
        dest_files.append(dest_file)
        dest_file_iter_bounds[dest_file] = [iter_0,iter_1]

    source_file_read_dict = {}
    source_file_read_index_sets = {}

    # # make dicts of the index and corresponding iters which will be read from these files
    for dest_file in dest_files:
        if print_level>=5:
            print('               - Creating reference dictionary for '+dest_file)

        source_file_read_dict[dest_file] = []
        source_file_read_index_sets[dest_file] = []

        dest_file_iter_0 = dest_file_iter_bounds[dest_file][0]
        dest_file_iter_1 = dest_file_iter_bounds[dest_file][1]
        if print_level>=5:
            print('                   - This file will cover iterations with start points '+str(dest_file_iter_0)+' to '+str(dest_file_iter_1))

        dest_file_iters = np.arange(dest_file_iter_0,dest_file_iter_1,iters_per_output)
        dest_file_iters += iters_per_output / 2
        if print_level >= 5:
            print('                   - This file will cover iterations with midpoints: '+str(np.min(dest_file_iters))+' to '+str(np.max(dest_file_iters)))

        source_file_names = []

        # get a list of files to use for each iteration
        for dest_iter in dest_file_iters:
            file_to_use = ''
            for file_name in file_names:
                iter_midpoints = iter_midpoint_dict[file_name]
                if dest_iter >= np.min(iter_midpoints) and dest_iter <= np.max(iter_midpoints):
                    if file_to_use == '':
                        file_to_use = file_name
                    else:
                        if int(file_to_use.split('.')[-2]) > int(file_name.split('.')[-2]):
                            file_to_use = file_name
            if file_to_use!='':
                source_file_names.append(file_to_use)

        if len(source_file_names)>0:

            # use the list of source file names to make a dict that lists which iters to read from it
            unique_list = list(set(source_file_names))
            if print_level >= 5:
                print('                    - Reading from file(s) '+', '.join(unique_list))
            for file_name in unique_list:
                source_file_read_dict[dest_file].append(file_name)
            source_file_read_dict[dest_file] = sorted(source_file_read_dict[dest_file])

            # loop through the source file names to identify which indices will be read from each one
            for file_name in source_file_read_dict[dest_file]:
                min_dest_iter_index = source_file_names.index(file_name)
                max_dest_iter_index = len(source_file_names) - 1 - source_file_names[::-1].index(file_name)
                min_dest_iter = dest_file_iters[min_dest_iter_index]
                max_dest_iter = dest_file_iters[max_dest_iter_index]
                # print('        - The file ' + file_name+ ' will cover iters '+str(min_dest_iter)+' to '+str(max_dest_iter))
                source_file_iters = iter_midpoint_dict[file_name]
                print(source_file_iters)
                print(min_dest_iter,max_dest_iter)
                # print('            - This file has iters '+str(np.min(source_file_iters))+' to '+str(np.max(source_file_iters)))
                # try:
                index_0 = list(source_file_iters).index(min_dest_iter)
                index_1 = list(source_file_iters).index(max_dest_iter)+1
                if file_name == source_file_read_dict[dest_file][-1]:
                    index_1 -= 1
                # print('            - This corresponds to indices '+str(index_0)+' ('+str(source_file_iters[index_0]) +') through '+str(index_1-1)+' ('+str(source_file_iters[index_1-1]) +')')
                source_file_read_index_sets[dest_file].append([index_0, index_1])
                # except:
                #     a=1
                #     # print('            - Didnt find the correct indices in this file')


    return(dest_files, source_file_read_dict, source_file_read_index_sets)

########################################################################################################################


def create_L2_BC_ref_file(config_dir, L2_model_name, parent_model_level, parent_model_name, print_level):

    if print_level>1:
        print('    - Creating the BC field reference for the '+L2_model_name)

    var_name = 'THETA'
    boundary = 'north'

    averaging_period = 3600
    if parent_model_level=='L1':
        seconds_per_iter = 300
    if parent_model_level=='L2':
        seconds_per_iter = 150

    # first read how many points are expected in each iter (read from the mask reference)
    mask_ref_file = os.path.join(config_dir, parent_model_level, parent_model_name,'input',
                                 'L1_dv_mask_reference_dict.nc')
    ds = nc4.Dataset(mask_ref_file)
    points_per_output = len(ds.groups['L3_'+boundary].variables['source_rows'][:])
    ds.close()

    Nr = 50

    # calculate the file name iterations to read from
    #     along with the iterations these files cover
    file_names, iter_midpoint_dict = find_dv_files_to_read(config_dir, parent_model_level, parent_model_name,
                                                           boundary, var_name, averaging_period, seconds_per_iter, points_per_output, Nr)
    if print_level>=2:
        print('    - Source file summary:')
    output = '{\n'
    for file_name in file_names:
        print('        - The file ' + file_name + ' has iterations centered on ' + str(np.min(iter_midpoint_dict[file_name])) +
              ' to ' + str(np.max(iter_midpoint_dict[file_name])))
        output += ' \''+file_name.split('.')[-2]+'\': '+'['+str(np.min(iter_midpoint_dict[file_name]))+\
                  ', '+str(np.max(iter_midpoint_dict[file_name]))+'],\n'
    output += '}'

    output_file = os.path.join(config_dir,'L2',L2_model_name,'input','L2_BC_source_ref.txt')
    f = open(output_file,'w')
    f.write(output)
    f.close()

    # calculate the destination file iteration bounds
    #     along with the source files data they will contain, and which indices of those source files they will contain
    dest_files, source_file_read_dict, source_file_read_index_sets = create_destination_file_list(config_dir, var_name, file_names,
                                                                                                  iter_midpoint_dict, averaging_period, seconds_per_iter, print_level)

    if print_level >= 2:
        print('    - Destination file summary:')
    output = '{\n'
    for file_name in dest_files:
        if print_level >= 2:
            print('        - The file ' + file_name + ' will be created from the following data:')
        output+=' \''+file_name.split('.')[-2]+'\': ['
        source_files = source_file_read_dict[file_name]
        index_sets = source_file_read_index_sets[file_name]
        iter_count = 0
        add_line = ''
        for s in range(len(source_files)):
            if len(index_sets)>0:
                if print_level >= 2:
                    print('            - From ' + source_files[s] + ', will read indices ' + str(
                    index_sets[s][0]) + ' though ' + str(index_sets[s][1]))
                iter_count += index_sets[s][1] - index_sets[s][0] + 1

                add_line += '[\'' + source_files[s].split('.')[-2] + '\', ' + '[' + str(index_sets[s][0]) + ', ' + str(index_sets[s][1]) + ']], '
        if int(iter_count / 24) == 1:
            output += add_line[:-2]
        else:
            output += ''
        if print_level >= 2:
            print('        - Total iterations for this file: '+str(iter_count)+' (= '+str(iter_count/24)+' days)')
        output+='],\n'
    output += '}'

    output_file = os.path.join(config_dir,parent_model_level, 'L2',L2_model_name,'input','L2_BC_dest_ref.txt')
    f = open(output_file,'w')
    f.write(output)
    f.close()



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1`, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_BC_fields_reference_dict(config_dir)