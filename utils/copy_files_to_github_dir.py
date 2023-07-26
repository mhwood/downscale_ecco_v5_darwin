

import argparse
import shutil
import os


def copy_files(config_dir, github_dir):

    level_names = ['L2','L1','L0']

    for level_name in level_names:
        if level_name not in os.listdir(github_dir):
            os.mkdir(os.path.join(github_dir,level_name))

    level_name_model_dict = {'L2':['L2_Disko_Bay', 'L2_Santa_Barbara', 'L2_Monterey_Bay'],
                             'L1':['L1_East_Pacific','L1_GOM','L1_W_Greenland','L1_mac_delta'],
                             'L0':[]}
    subdirs = []
    utils_subdirs = ['init_file_creation','plot_creation']

    for level_name in level_names:

        print('    - Copying the models in level '+level_name)

        #################################################################
        # copy the model files first
        model_names = level_name_model_dict[level_name]

        for model_name in model_names:

            print('        - Copying the files for model ' + model_name)

            if model_name not in os.listdir(os.path.join(github_dir,level_name)):
                os.mkdir(os.path.join(github_dir,level_name,model_name))

            #####################################
            # copy the code and namelist files
            for dir_name in subdirs:
                if dir_name in os.listdir(os.path.join(config_dir,level_name,model_name)):
                    if dir_name not in os.listdir(os.path.join(github_dir,level_name,model_name)):
                        os.mkdir(os.path.join(github_dir,level_name,model_name,dir_name))
                    for file_name in os.listdir(os.path.join(config_dir,level_name,model_name,dir_name)):
                        if file_name[0]!='.' and file_name[-3:]!='.nc':
                            shutil.copyfile(os.path.join(config_dir,level_name,model_name,dir_name,file_name),
                                            os.path.join(github_dir,level_name,model_name,dir_name,file_name))

            #####################################
            # copy the utils files
            if 'utils' not in os.listdir(os.path.join(github_dir, level_name, model_name)):
                os.mkdir(os.path.join(github_dir, level_name, model_name, 'utils'))
            for dir_name in utils_subdirs:
                if dir_name in os.listdir(os.path.join(config_dir, level_name, model_name, 'utils')):
                    if dir_name not in os.listdir(os.path.join(github_dir, level_name, model_name,'utils')):
                        os.mkdir(os.path.join(github_dir, level_name, model_name, 'utils', dir_name))
                    for file_name in os.listdir(os.path.join(config_dir, level_name, model_name, 'utils',dir_name)):
                        if file_name[0] != '.' and file_name[-3:]=='.py' and file_name[0] != '_' and file_name!='init_files' and file_name!='results':
                            shutil.copyfile(os.path.join(config_dir, level_name, model_name, 'utils',dir_name, file_name),
                                            os.path.join(github_dir, level_name, model_name, 'utils',dir_name, file_name))
                    if dir_name == 'plot_creation' and 'init_files' in os.listdir(os.path.join(config_dir, level_name, model_name,'utils', dir_name)):
                        if 'init_files' not in os.listdir(os.path.join(github_dir, level_name, model_name,'utils',dir_name)):
                            os.mkdir(os.path.join(github_dir, level_name, model_name, 'utils', dir_name,'init_files'))
                        for file_name in os.listdir(os.path.join(config_dir, level_name, model_name, 'utils',dir_name,'init_files')):
                            if file_name[0] != '.' and file_name[-3:]=='.py' and file_name[0] != '_':
                                shutil.copyfile(os.path.join(config_dir, level_name, model_name, 'utils',dir_name,'init_files', file_name),
                                                os.path.join(github_dir, level_name, model_name, 'utils',dir_name,'init_files', file_name))
                    if dir_name == 'plot_creation' and 'results' in os.listdir(os.path.join(config_dir, level_name, model_name,'utils', dir_name)):
                        if 'results' not in os.listdir(os.path.join(github_dir, level_name, model_name,'utils',dir_name)):
                            os.mkdir(os.path.join(github_dir, level_name, model_name, 'utils', dir_name,'results'))
                        for file_name in os.listdir(os.path.join(config_dir, level_name, model_name, 'utils',dir_name,'results')):
                            if file_name[0] != '.' and file_name[-3:]=='.py' and file_name[0] != '_':
                                shutil.copyfile(os.path.join(config_dir, level_name, model_name, 'utils',dir_name,'results', file_name),
                                                os.path.join(github_dir, level_name, model_name, 'utils',dir_name,'results', file_name))
            for file_name in os.listdir(os.path.join(config_dir, level_name, model_name, 'utils')):
                if file_name[-3:]=='.sh' or file_name[-3:]=='.py':
                    shutil.copyfile(os.path.join(config_dir, level_name, model_name, 'utils', file_name),
                                    os.path.join(github_dir, level_name, model_name, 'utils', file_name))

        #####################################
        # copy the general level utils files
        print('        - Copying the utils for level ' + level_name)
        if 'utils' in os.listdir(os.path.join(config_dir, level_name)):
            if 'utils' not in os.listdir(os.path.join(github_dir, level_name)):
                os.mkdir(os.path.join(github_dir, level_name, 'utils'))
            for dir_name in utils_subdirs:
                if dir_name in os.listdir(os.path.join(config_dir, level_name, 'utils')):
                    if dir_name not in os.listdir(os.path.join(github_dir, level_name, 'utils')):
                        os.mkdir(os.path.join(github_dir, level_name, 'utils', dir_name))
                    for file_name in os.listdir(os.path.join(config_dir, level_name, 'utils', dir_name)):
                        if file_name[0] != '.' and file_name[0] != '_' and file_name!='init_files' and file_name!='results':
                            shutil.copyfile(
                                os.path.join(config_dir, level_name, 'utils', dir_name, file_name),
                                os.path.join(github_dir, level_name, 'utils', dir_name, file_name))
                    if dir_name == 'plot_creation' and 'init_files' in os.listdir(os.path.join(config_dir, level_name, 'utils',dir_name)):
                        if 'init_files' not in os.listdir(os.path.join(github_dir, level_name, 'utils',dir_name)):
                            os.mkdir(os.path.join(github_dir, level_name, 'utils', dir_name, 'init_files'))
                        for file_name in os.listdir(os.path.join(config_dir, level_name, 'utils', dir_name,'init_files')):
                            if file_name[0] != '.' and file_name[0] != '_':
                                shutil.copyfile(
                                    os.path.join(config_dir, level_name, 'utils', dir_name,'init_files', file_name),
                                    os.path.join(github_dir, level_name, 'utils', dir_name,'init_files', file_name))
                    if dir_name == 'plot_creation' and 'results' in os.listdir(os.path.join(config_dir, level_name, 'utils',dir_name)):
                        if 'results' not in os.listdir(os.path.join(github_dir, level_name, 'utils',dir_name)):
                            os.mkdir(os.path.join(github_dir, level_name, 'utils', dir_name, 'results'))
                        for file_name in os.listdir(os.path.join(config_dir, level_name, 'utils', dir_name,'results')):
                            if file_name[0] != '.' and file_name[0] != '_':
                                shutil.copyfile(
                                    os.path.join(config_dir, level_name, 'utils', dir_name,'results', file_name),
                                    os.path.join(github_dir, level_name, 'utils', dir_name,'results', file_name))

            for file_name in os.listdir(os.path.join(config_dir, level_name, 'utils')):
                if file_name[-3:] == '.sh' or file_name[-3:] == '.py':
                    shutil.copyfile(os.path.join(config_dir, level_name, 'utils', file_name),
                                    os.path.join(github_dir, level_name, 'utils', file_name))

    ####################################################################################################################
    # copy the general code and namelist files
    print('    - Copying the general utils ')
    if 'utils' not in os.listdir(github_dir):
        os.mkdir(os.path.join(github_dir, 'utils'))
    for dir_name in utils_subdirs:
        if dir_name in os.listdir(os.path.join(config_dir, 'utils')):
            if dir_name not in os.listdir(os.path.join(github_dir, 'utils')):
                os.mkdir(os.path.join(github_dir, 'utils', dir_name))
            for file_name in os.listdir(os.path.join(config_dir, 'utils', dir_name)):
                if file_name[0] != '.' and file_name[0] != '_':
                    shutil.copyfile(
                        os.path.join(config_dir,  'utils', dir_name, file_name),
                        os.path.join(github_dir,  'utils', dir_name, file_name))
    for file_name in os.listdir(os.path.join(config_dir,  'utils')):
        if file_name[-3:] == '.sh' or file_name[-3:] == '.py' and file_name[0] != '_' and file_name!='copy_files_to_github_dir.py':
            shutil.copyfile(os.path.join(config_dir, 'utils', file_name),
                            os.path.join(github_dir, 'utils', file_name))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-g", "--github_dir", action="store",
                        help="The directory which will be pushed to Github.", dest="github_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    github_dir = args.github_dir

    copy_files(config_dir, github_dir)
