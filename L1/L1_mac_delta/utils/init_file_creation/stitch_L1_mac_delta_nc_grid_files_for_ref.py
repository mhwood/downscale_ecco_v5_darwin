
import os
import numpy as np
import netCDF4 as nc4
import argparse
import ast

########################################################################################################################

def read_L1_grid_tile_geometry(config_dir, model_name, Nr, sNx, sNy, ordered_nonblank_tiles, ordered_tile_rotations):

    stitched_XC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_YC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    stitched_XG = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_YG = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    stitched_AngleCS = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_AngleSN = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    stitched_DXC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0]) + 1))
    stitched_DYC = np.zeros((sNy * len(ordered_nonblank_tiles) + 1, sNx * len(ordered_nonblank_tiles[0])))

    stitched_Depth = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_rA = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    stitched_hFacC = np.zeros((Nr, sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_hFacS = np.zeros((Nr, sNy * len(ordered_nonblank_tiles) + 1, sNx * len(ordered_nonblank_tiles[0])))
    stitched_hFacW = np.zeros((Nr, sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0]) + 1))

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            rotations = ordered_tile_rotations[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    XG = ds.variables['XG'][:, :]
                    YG = ds.variables['YG'][:, :]
                    # AngleCS = ds.variables['AngleCS'][:, :]
                    # AngleSN = ds.variables['AngleSN'][:, :]
                    DXC = ds.variables['dxC'][:, :]
                    DYC = ds.variables['dyC'][:, :]
                    hFacC = ds.variables['HFacC'][:, :, :]
                    hFacS = ds.variables['HFacS'][:, :, :]
                    hFacW = ds.variables['HFacW'][:, :, :]
                    DRF = ds.variables['drF'][:]
                    Depth = ds.variables['Depth'][:, :]
                    rA = ds.variables['rA'][:, :]
                    ds.close()

                    XG = XG[:-1,:-1]
                    YG = YG[:-1,:-1]

                    for i in range(rotations):
                        XC = np.rot90(XC)
                        YC = np.rot90(YC)
                        XG = np.rot90(XG)
                        YG = np.rot90(YG)
                        # AngleCS = np.rot90(AngleCS)
                        # AngleSN = np.rot90(AngleSN)

                        Depth = np.rot90(Depth)
                        rA = np.rot90(rA)

                        DXC = np.rot90(DXC)
                        DYC = np.rot90(DYC)

                        hFacC = np.rot90(hFacC, axes=(1, 2))
                        hFacS = np.rot90(hFacS, axes=(1, 2))
                        hFacW = np.rot90(hFacW, axes=(1, 2))

                    stitched_XC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = XC
                    stitched_YC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = YC

                    stitched_XG[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = XG
                    stitched_YG[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = YG

                    # stitched_AngleCS[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = AngleCS
                    # stitched_AngleSN[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = AngleSN

                    stitched_DXC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx+1] = DXC
                    stitched_DYC[r * sNy:(r + 1) * sNy+1, c * sNx:(c + 1) * sNx] = DYC

                    stitched_hFacC[:, r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = hFacC
                    stitched_hFacS[:, r * sNy:(r + 1) * sNy + 1, c * sNx:(c + 1) * sNx] = hFacS
                    stitched_hFacW[:, r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx + 1] = hFacW

                    stitched_Depth[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = Depth
                    stitched_rA[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = rA

    return(stitched_XC, stitched_YC, stitched_XG, stitched_YG, stitched_AngleCS, stitched_AngleSN, stitched_DXC, stitched_DYC,
           stitched_hFacC, stitched_hFacS, stitched_hFacW,
           stitched_Depth, stitched_rA, DRF)


def write_grid_to_nc(config_dir, model_name,
                     XC, YC, XG, YG, AngleCS, AngleSN, DXC, DYC, hFacC, hFacS, hFacW, DRF, Depth, rA):

    output_path = os.path.join(config_dir, 'nc_grids', model_name+'_grid.nc')

    ds = nc4.Dataset(output_path,'w')

    ds.createDimension('X', np.shape(XC)[1])
    ds.createDimension('Y', np.shape(XC)[0])
    ds.createDimension('Xp1', np.shape(DXC)[1])
    ds.createDimension('Yp1', np.shape(DYC)[0])
    ds.createDimension('Z',np.shape(DRF)[0])

    var = ds.createVariable('XC','f4',('Y','X'))
    var[:,:] = XC

    var = ds.createVariable('YC', 'f4', ('Y', 'X'))
    var[:, :] = YC

    var = ds.createVariable('XG', 'f4', ('Y', 'X'))
    var[:, :] = XG

    var = ds.createVariable('YG', 'f4', ('Y', 'X'))
    var[:, :] = YG

    var = ds.createVariable('AngleCS', 'f4', ('Y', 'X'))
    var[:, :] = AngleCS

    var = ds.createVariable('AngleSN', 'f4', ('Y', 'X'))
    var[:, :] = AngleSN

    var = ds.createVariable('dxC', 'f4', ('Y', 'Xp1'))
    var[:, :] = DXC

    var = ds.createVariable('dyC', 'f4', ('Yp1', 'X'))
    var[:, :] = DYC

    var = ds.createVariable('HFacC', 'f4', ('Z', 'Y', 'X'))
    var[:, :, :] = hFacC

    var = ds.createVariable('HFacW', 'f4', ('Z', 'Y', 'Xp1'))
    var[:, :, :] = hFacW

    var = ds.createVariable('HFacS', 'f4', ('Z', 'Yp1', 'X'))
    var[:, :, :] = hFacS

    var = ds.createVariable('drF', 'f4', ('Z',))
    var[:] = DRF

    var = ds.createVariable('Depth', 'f4', ('Y', 'X'))
    var[:, :] = Depth

    var = ds.createVariable('rA', 'f4', ('Y', 'X'))
    var[:, :] = rA

    ds.close()

def stitch_grid_files(config_dir):

    print('Stitching the nc grid files')

    model_name = 'L1_mac_delta'
    ordered_tiles = [[1,2,3,4,5,6],[7,8,9,10,11,12],[13,14,15,16,17,18],[19,20,21,22,23,24]]
    ordered_tile_rotations = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
    Nr = 80
    sNx = 204
    sNy = 186

    XC, YC, XG, YG, AngleCS, AngleSN, DXC, DYC, hFacC, hFacS, hFacW, Depth, rA, DRF \
        = read_L1_grid_tile_geometry(config_dir, model_name, Nr, sNx, sNy, ordered_tiles, ordered_tile_rotations)

    # print('   - Copying interior AngleCS and AngleSN values to boundary because these values are messed up
    AngleCS[-1, :] = AngleCS[-2, :]
    AngleCS[0, :] = AngleCS[1, :]
    AngleCS[:, -1] = AngleCS[:, -2]
    AngleSN[-1, :] = AngleSN[-2, :]
    AngleSN[0, :] = AngleSN[1, :]
    AngleSN[:, -1] = AngleSN[:, -2]

    # print('   - Copying interior HFacS values to boundary because these values are messed up
    hFacS[:, -1, :] = hFacS[:, -2, :]
    hFacS[:, 0, :] = hFacS[:, 1, :]

    # print('   - Copying interior HFacS values to boundary because these values are messed up
    hFacW[:, :, -1] = hFacW[:, :, -2]
    hFacW[:, :, 0] = hFacW[:, :, 1]

    write_grid_to_nc(config_dir, model_name,
                     XC, YC, XG, YG, AngleCS, AngleSN, DXC, DYC, hFacC, hFacS, hFacW, DRF, Depth, rA)

    # zero_rows = 0
    # for i in range(np.shape(hFacC)[0]):
    #     if np.all(hFacC[i,:,:]==0):
    #         zero_rows+=1
    #
    # if zero_rows>1:
    #     print('    - The grid has '+str(zero_rows)+' zero rows - consider chopping off some of them to reduce computational time')




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    stitch_grid_files(config_dir)