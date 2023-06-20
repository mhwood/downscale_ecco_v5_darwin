import os
import numpy as np
import simplegrid as sg
import shapefile
import argparse


def read_ecco_geometry_to_faces(ecco_dir,llc):
    # grid_keys = ['XC', 'YC', 'DXF', 'DYF', 'RAC', 'XG', 'YG', 'DXV', 'DYU',
    #              'RAZ', 'DXC', 'DYC', 'RAW', 'RAS', 'DXG', 'DYG']
    grid_file_dir = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles')
    XC_faces = {}
    YC_faces = {}

    XG_faces = {}
    YG_faces = {}
    for i in [1, 2, 3, 4, 5]:
        if i < 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc,
                                                   3 * llc)
        if i == 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc,
                                                   llc)
        if i > 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), 3 * llc,
                                                   llc)
        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face

        XG_face = grid_dict['XG'].T
        YG_face = grid_dict['YG'].T
        XG_faces[i] = XG_face
        YG_faces[i] = YG_face

    return(XC_faces,YC_faces,XG_faces,YG_faces)

def read_ecco_depth_to_faces(ecco_dir,llc):

    file_path = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init','bathy_llc270')

    grid_array = np.fromfile(file_path,'>f4')

    field_faces = {}

    points_counted = 0
    for i in range(1, 6):
        if i < 3:
            n_points = 3 * llc * llc
            grid = grid_array[points_counted:points_counted + n_points]
            grid = np.reshape(grid, (3 * llc, llc))
        if i == 3:
            n_points = llc * llc
            grid = grid_array[points_counted:points_counted + n_points]
            grid = np.reshape(grid, (llc, llc))
        if i > 3:
            n_points = 3 * llc * llc
            grid = grid_array[points_counted:points_counted + n_points]
            grid = np.reshape(grid, (llc, 3 * llc))
        field_faces[i] = grid
        points_counted += n_points

    return(field_faces)

def face_grid_to_shp(output_file,face,XC,YC,XG,YG,depth):
    sf = shapefile.Writer(output_file)

    sf.field('face', 'N')
    sf.field('row', 'N')
    sf.field('col', 'N')
    sf.field('XC', 'F',15,15)
    sf.field('YC', 'F',15,15)
    sf.field('Depth', 'N')

    if face==3:
        for row in range(np.shape(XG)[0] - 1):
            for col in range(np.shape(XG)[1] - 1):
                if int(XG[row,col]) < 0:
                    polygon = [[XG[row, col], YG[row, col]],
                               [XG[row + 1, col], YG[row + 1, col]],
                               [XG[row + 1, col + 1], YG[row + 1, col + 1]],
                               [XG[row, col + 1], YG[row, col + 1]]]
                    sf.poly([polygon])
                    print(face,row, col,XC[row, col], YC[row, col])
                    sf.record(face,row, col,XC[row, col], YC[row, col],np.round(depth[row,col]))
    else:
        for row in range(np.shape(XG)[0] - 1):
            for col in range(np.shape(XG)[1] - 1):
                # if int(Depth[row, col]) > 0:
                # polygon = [[XG[row,col]-XG[row,col],YG[row,col]-YG[row,col]],
                #            [XG[row,col]+XG[row,col],YG[row,col]-YG[row,col]],
                #            [XG[row,col]+XG[row,col],YG[row,col]+YG[row,col]],
                #            [XG[row,col]-XG[row,col],YG[row,col]+YG[row,col]]]
                polygon = [[XG[row, col], YG[row, col]],
                           [XG[row + 1, col], YG[row + 1, col]],
                           [XG[row + 1, col + 1], YG[row + 1, col + 1]],
                           [XG[row, col + 1], YG[row, col + 1]]]
                sf.poly([polygon])
                sf.record(face,row, col,XC[row, col], YC[row, col],np.round(depth[row,col]))

    sf.close()

    f = open(output_file+'.prj','w')
    f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]')
    f.close()
    a=1

########################################################################################################################

def tiles_to_shp(config_dir,ecco_path):
    llc = 270

    XC_faces, YC_faces, XG_faces, YG_faces = read_ecco_geometry_to_faces(ecco_path,llc)

    depth_faces = read_ecco_depth_to_faces(ecco_path,llc)

    if 'grid_shapefiles' not in os.listdir(os.path.join(config_dir,'L0','input')):
        os.mkdir(os.path.join(config_dir,'L0','input','grid_shapefiles'))

    for face in [1,2,3,4,5]:
        output_file = os.path.join(config_dir,'L0','input','grid_shapefiles','face_'+str(face)+'_grid')
        face_grid_to_shp(output_file, face, XC_faces[face], YC_faces[face], XG_faces[face], YG_faces[face], depth_faces[face])







if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_path = args.ecco_path

    tiles_to_shp(config_dir,ecco_path)
