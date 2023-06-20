

import numpy as np


def calculate_flux_into_domain(boundary, var_name,
                                   velocity, tracer,
                                   dxG, dyG, drF, HFacS, HFacW,
                                   scale_factor = 1, density = 1024, apply_density = True):

    heat_capacity = 3974  # J/kg/C

    if boundary=='north':
        HFac = HFacS[:,-2,:]
        width = dxG[-2,:]
    if boundary=='south':
        HFac = HFacS[:,1,:]
        width = dxG[1, :]
    if boundary=='east':
        HFac = HFacW[:,:,-2]
        width = dyG[:, -2]
    if boundary=='west':
        HFac = HFacW[:,:,1]
        width = dyG[:, 1]

    flux_grid = np.zeros((np.shape(velocity)[0],len(drF),len(width)))

    if not apply_density:
        density = 1

    for time_step in range(np.shape(velocity)[0]):
        for k in range(len(drF)):
            if var_name in ['UVEL','VVEL']:
                flux_grid[time_step, k, :] += width * velocity[time_step, k, :] * drF[k] * HFac[k, :] * scale_factor
            elif var_name in ['Theta']:
                flux_grid[time_step, k, :] += width * velocity[time_step, k, :] * drF[k] * HFac[k,:] *\
                                              heat_capacity * density * tracer[time_step,k,:] * scale_factor
            else :
                flux_grid[time_step, k, :] += width * velocity[time_step, k, :] * drF[k] * HFac[k,:] *\
                                              density * tracer[time_step,k,:] * scale_factor

    if boundary=='north':
        flux_grid *= -1
    if boundary=='east':
        flux_grid *= -1

    return(flux_grid)

