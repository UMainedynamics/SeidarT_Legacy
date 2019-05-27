#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created by Steven Bernsen
We are given a range of velocities of different materials found empirically. 
For isotropic materials we can determine the Lame constants from the equations:
    Vp = sqrt( lambda + 2 mu  / rho ),
    Vs = sqrt( mu / rho ),
    c11 = c22 = c33 = lambda + 2 mu,
    c12 = c13 = c23 = lambda,
    c44 = c55 = c66 = mu
"""

import numpy as np


# =============================================================================
#                       Define material dictionaries
# =============================================================================

"""
Seismic values can be found in: 
      Acoustics of Porous Media (1992), Bourbie, Coussy, and Zinszner
      
Permittivity values can be found in:
        Electrical Properties of Rocks and Minerals
    The following values are provided by:
        Seismic Velocities - https://pangea.stanford.edu/courses/gp262/Notes/8.SeismicVelocity.pdf
        Permeabilitues - http://www.geo.umass.edu/faculty/wclement/dielec.html
        Conductivities - Duba et al. (1977), Duba et al. (1978), Watanabe (1970), 
                    Mohammadi and Mohammadi (2016),
                    https://www.nrcs.usda.gov/INTERNET/FSE_DOCUMENTS/NRCS142P2_053280.PDF,
                    https://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity
    Values are: 
        Vp_min, Vp_max, Vs_min, Vs_max, Rel_Perm_min, Rel_Perm_max, Conductivity_min, Conductivity_max
    
    Permittivity is given as the relative permittivity and for most cases we 
    will assume that relative permeability is unity; however, if we include
    materials that are high in magnetite, hematite, etc. then we will need to
    accomadate for better permeability estimates.
    """


isotropic_materials = {
    "air":np.array([343, 343, 0.0, 0.0, 1.0, 1.0, 1.0e-16, 1.0e-15]),
    "ice1h":np.array([3400, 3800, 1700, 1900, 3.1, 3.22, 1.0e-7, 1.0e-6]),
    "soil":np.array([300, 700, 100, 300, 3.9, 29.4, 1.0e-2, 1.0e-1]),
    "water":np.array([1450, 1500, 0, 0, 80.36, 80.36, 5.5e-6, 5.0e-2]), # This can change drastically depending on the ions in solution
    "oil":np.array([1200, 1250, 0, 0, 2.07, 2.14, 5.7e-8, 2.1e-7]),
    "dry_sand":np.array([400, 1200, 100, 500, 2.9, 4.7, 1.0e-3, 1.0e-3]), # perm porositiy dependence
    "wet_sand":np.array([1500, 2000, 400, 600, 2.9, 105, 2.5e-4, 1.2e-3]), 
    "granite":np.array([4500, 6000, 2500, 3300, 4.8, 18.9, 4.0e-5, 2.5e-4]),
    "gneiss":np.array([4400, 5200, 2700, 3200, 8.5, 8.5, 2.5e-4, 2.5e-3]),
    "basalt":np.array([5000, 6000, 2800, 3400, 12, 12, 1.0e-6, 1.0e-4]),
    "limestone":np.array([3500, 6000, 2000, 3300, 7.8, 8.5, 2.5e-4, 1.0e-3]),
    "anhydrite":np.array([4000, 5500, 2200, 3100, 5, 11.5, 1.0e-6, 1.0e-5]), # permittivity value from Gypsum
    "coal":np.array([2200, 2700, 1000, 1400, 5.6, 6.3, 1.0e-8, 1.0e-3]), # This has a water dependency
    "salt":np.array([4500, 5500, 2500, 3100, 5.6, 5.6, 1.0e-7, 1.0e2]) # This is dependent on ions and water content
}


# ======================================
# Just a check to make sure it's working
# def print_vals(self):
#     pass
    # p = ice_permittivity(-2)
    # s = ice_stiffness(-10, 3)
    # print(s)
# ======================================


# -----------------------------------------------------------------------------
def pressure_array(im, temp, rho, dz, pore = 0, wc = 0):
    # im is an m-by-n array of integer values corresponding to the material
    # temp and rho are k-by-1 indices


    # First match the size of 
    k = np.unique(im)

    # if len(pore) > 1 and not len(pore) == k :
    #     pore = np.ones(k)*pore

    # if len(wc) > 1 and not len(wc) == k:
    #     wc = np.ones(k)*pore


    m, n = im.shape
    # allocate the pressure, temperature and density
    pressure = np.zeros([m, n])
    density = np.zeros([m, n])

    if not temp.shape == im.shape:
        temperature = np.zeros([m, n])


    for j in range(0,n):
        for i in range(0, m):
            temperature[i,j] = temp[ im[i,j] ]
            density[i,j] = porewater_correction(temperature[i,j], rho[ im[i,j] ], 
                pore[ im[i,j]], wc[ im[i,j]])
            pressure[i,j] = density[i,j] * 9.80665 * i * dz
            
    
    return(temperature, density, pressure)


def anisotrpic_boolean(im, matbool, angvect):

    m,n = im.shape

    abool = np.zeros([m,n], dtype = bool)
    afile = np.zeros([m,n], dtype = str)

    for i in range(0, m):
        for j in range(0,n):

            if matbool[ im[i,j] ] == 'True':
                abool[i,j] = True
                afile[i,j] = angvect[ im[i,j] ]

    return(abool, afile)

# =============================================================================
# -----------------------------------------------------------------------------
def get_seismic(material_name = None, temp=None, rho=None, pore = 0, 
    wc = 0, abool=None, angfile=None):

    m = len(temp)
    tensor = np.zeros([m, 11])

    # Adjust the stiffness tensor according to the pressure and temperature conditions for ice
    for ind in range(0, m):
        density = porewater_correction(temp[ind], rho[ind], pore[ind], wc[ind] )

        if abool[ind] and material_name[ind] == 'ice1h':
            euler = read_ang(angfile[ind])
            p = len(euler[:,0])

            cvoigt = np.zeros([6,6])
            creuss = np.zeros([6,6])
            C = np.zeros([6,6])
            
            # Assume a constant pressure of 0.1 MPa
            C = ice_stiffness(temp[ind], 0.1)
            S = np.linalg.inv(C)

            for k in range(0, p):
                R = rotator_zxz(euler[k,:] )
                M = bond(R)
                N = np.linalg.inv(M)
                cvoigt = cvoigt + ( np.matmul( M, np.matmul(C, M.T) ) )
                creuss = creuss + ( np.matmul( N, np.matmul(S, N.T) ) )

            cvoigt = cvoigt/p
            creuss = creuss/p 
            creuss = np.linalg.inv(creuss) 
            
            # Calculate the hill average 
            C = (cvoigt + creuss)/2

        else:
            material_limits = isotropic_materials[ material_name[ind] ]
            C = isotropic_stiffness_tensor(0.1, density, material_limits )

        tensor[ind, :] = ( ind, C[0,0], C[0,1], C[0,2], C[1,1], C[1,2], C[2,2], C[3,3], C[4,4], C[5,5], density)

    return(tensor)


# -----------------------------------------------------------------------------
def get_perm(material_name = None, temp=None, pore = 0, 
    wc = 0, abool=None, angfile=None):

    m = len(temp)
    tensor = np.zeros([m, 7])

    # Adjust the stiffness tensor according to the pressure and temperature conditions for ice
    for ind in range(0, m):
        
        cond = isotropic_permittivity_tensor(temp[ind], pore[ind], wc[ind], material_name[ind])[1]
        if abool[ind] and material_name[ind] == 'ice1h':
            euler = read_ang(angfile[ind])
            p = len(euler[:,0])

            pvoigt = np.zeros([3,3])
            preuss = np.zeros([3,3])
            P = np.zeros([3,3])
            
            # Assume a constant pressure of 0.1 MPa
            P = ice_permittivity(temp[ind])
            S = np.linalg.inv(P)

            for k in range(0, p):
                R = rotator_zxz(euler[k,:] )
            
                Ri = np.linalg.inv(R)
                pvoigt = pvoigt + ( np.matmul( R, np.matmul(P, R.T) ) )
                preuss = preuss + ( np.matmul( Ri, np.matmul(S, Ri.T) ) )

            pvoigt = pvoigt/p
            preuss = preuss/p 
            preuss = np.linalg.inv(preuss) 
            
            # Calculate the hill average 
            P = (pvoigt + preuss)/2

        else:
            P = np.round(isotropic_permittivity_tensor(temp[ind], pore[ind], wc[ind], material_name[ind])[0], 3)

        tensor[ind, :] = np.array([ind, P[0,0], P[1,1], P[2,2], cond[0,0], cond[1,1], cond[2,2] ] )
        

    return(tensor)


# -----------------------------------------------------------------------------
# =============================================================================


def isotropic_stiffness_tensor(pressure, density, material_limits):

    Vp = material_limits[0:2]
    Vs = material_limits[2:4]
    cp = 2*(Vp[1] - Vp[0])/np.pi 
    cs = 2*(Vs[1] - Vs[0])/np.pi
    
    # Correct for pressure
    pvelocity = cp*np.arctan(pressure ) + Vp[0]
    svelocity = cs*np.arctan(pressure) + Vs[0]
    
    # Compute the lame parameters
    mu = density*(svelocity**2)
    lam = density*(pvelocity**2) - 2*mu

    # Assign the matrix
    C = np.zeros([6,6])
    C[0:3,0:3] = lam
    np.fill_diagonal(C, C.diagonal() + mu)
    C[0,0] = lam + 2*mu
    C[1,1]= C[1,1] + mu
    C[2,2] = C[2,2] + mu

    return(C)

# -----------------------------------------------------------------------------
def isotropic_permittivity_tensor(temperature, porosity, water_content, material_name):

    material_limits = isotropic_materials[ material_name ]
    perm0 = material_limits[4]
    perm1 = material_limits[5]

    cond0 = material_limits[6]
    cond1 = material_limits[7]

    # Calculate the slope         
    if material_name == 'ice1h':
        # We'll assume that the maximum porosity of ice (a.k.a. fresh pow pow)
        # is 85%. The porosity is given as percent [0,100]
        perm0 = 3.1884 + 9.1e-4 * temperature
        perm_coef = (perm1 - perm0)/85
        cond_coef = (cond1 - cond0)/85
        permittivity = np.eye(3,3) * (perm_coef*(porosity) + perm0)
        conductivity = np.eye(3,3) * (cond_coef*(porosity) + cond0)
            
    elif material_name == 'soil' or material_name == 'dry sand':
        # The limit of these two materials is around 55%
        perm_coef = (perm1 - perm0)/55
        cond_coef = (cond1 - cond0)/55
        permittivity = np.eye(3,3) * (perm_coef*(porosity) + perm0)
        conductivity = np.eye(3,3) * (cond_coef*(porosity) + cond0)
    
    elif material_name == 'salt':
        # Permittivity will change a little bit but let's neglect it for now
        permittivity = np.eye(3,3) * perm0
        # Let's use a simple linear trend for a maximum of 20%? water content
        cond_coef = (cond1 - cond0)/20
        conductivity = np.eye(3,3) * (cond_coef*(water_content) +  cond0 )
    
    elif material_name == 'water' or material_name == 'oil':
        # Water and oil do not have a porosity.
        permittivity = np.eye(3,3) * material_limits[4]
        conductivity = np.eye(3,3) * material_limits[6]
    else:
        # For other materials we'll assume the limit is around 3%
        perm_coef = (perm1 - perm0)/3
        cond_coef = (cond1 - cond0)/3
        permittivity = np.eye(3,3) * (perm_coef*(porosity) + perm0)
        conductivity = np.eye(3,3) * (cond_coef*(porosity) + cond0)
    
    return(permittivity, conductivity)

# -----------------------------------------------------------------------------
def porewater_correction(temperature, density, porosity, water_content ):

    rho_air = 0.02897/(8.2057338e-5 * (273 + temperature) )
    rho_water = -4.6074e-7*temperature**4 + 1.0326e-4*temperature**3 - 1.0833e-2*temperature**2 + 9.4207e-2*temperature**1 + 999.998

    # There are certain limits such as phase changes so let's put practical limits on this
    rho_water = np.max( (rho_water, 950) )
    # the water content is the percent of pores that contain water
    rho_wc = (1-water_content/100)*rho_air + (water_content/100)*rho_water
    density = (1-porosity/100)*density + (porosity/100)

    return(density)

# -----------------------------------------------------------------------------
def ice_stiffness(temperature = None, pressure = 0) :

    # Allocate space for the stiffness tensor
    C = np.zeros([6,6])
    
    C[0,0] = 136.813 - 0.28940*temperature - 0.00178270*(temperature**2) \
      + 4.6648*pressure - 0.13501*(pressure**2) 
    C[0,1] = 69.4200 - 0.14673*temperature - 0.00090362*(temperature**2) \
      + 5.0743*pressure + .085917*(pressure**2)
    C[0,2] = 56.3410 - 0.11916*temperature - 0.00073120*(temperature**2) \
      + 6.4189*pressure - .52490*(pressure**2)
    C[2,2] = 147.607 - 0.31129*temperature - 0.0018948*(temperature**2) \
      + 4.7546*pressure - .11307*(pressure**2)
    C[3,3] = 29.7260 - 0.062874*temperature - 0.00038956*(temperature**2) \
      + 0.5662*pressure + .036917*(pressure**2)
    
    # Fill in the symmetry
    C[1,1] = C[0,0]
    C[1,0] = C[0,1]
    C[2,0] = C[0,2]
    C[1,2] = C[0,2]
    C[2,1] = C[1,2]
    C[4,4] = C[3,3]
    C[5,5] = (C[0,0] - C[0,1] )/2
    
    stiffness = C*1e8

    return(stiffness)

# -----------------------------------------------------------------------------

def ice_permittivity(temperature):
    #Allocate 
    P = np.zeros([3,3])

    # The following is for 2-10 GHz. The details can be found in 
    # Fujita et al. (2000)
    perm = 3.1884 + 9.1e-4 * temperature
    dP = 0.0256 + 3.57e-5 * (6.0e-6) * temperature
    permittivity = np.eye(3,3) * perm 
    permittivity[2,2] = perm + dP 

    return(permittivity)


# -----------------------------------------------------------------------------

def read_ang(filepath):
    """
    The input .ang file will have the columns as 
        c1-3    Euler angles (radians; Bunge's notation - z-x-z rotation )
        c4,5    horizontal, vertical respectively
        c6      image quality
        c7      confidence index
        c8      phase ID
        c9      detector intensity
        c10     fit
    Refer to this thread for more description of the aforementioned
        https://www.researchgate.net/post/How_can_I_get_X_Y_position_data_of_a_single_grain_from_EBSD_scan
    """
    
    # Load the file in as a data frame
    euler = np.genfromtxt(filepath, delimiter = " ")

    # take only the euler angles...for now
    if euler.shape[0] > 3 :
        euler = euler[:,0:3]

    # Unfortunately, the space delimiters are inconsistent :(
    # We know there are 10 columns and the rows shouldn't contain all NA's
    m, n = np.shape(euler)

    # reshape to M-by-1 vector
    euler = euler.reshape(m*n,1)

    # remvoe 'nan'
    euler = euler[~np.isnan(euler)]

    # reshape back to array
    euler = euler.reshape(m, int( len(euler)/m ) )

    # save ferris
    return(euler)


# -----------------------------------------------------------------------------
def rotator_zxz(eul):
    # From the 3 euler angles for the zxz rotation, compute the rotation matrix
    R = np.zeros([3,3])
    D = np.zeros([3,3])
    C = np.zeros([3,3])
    B = np.zeros([3,3])

    D[0,:] = [ np.cos( eul[0] ), -np.sin( eul[0] ), 0.0 ]
    D[1,:] = [ np.sin( eul[0] ), np.cos( eul[0] ), 0.0 ]
    D[2,:] = [ 0.0, 0.0, 1.0 ]

    C[0,:] = [ 1.0, 0.0, 0.0 ]
    C[1,:] = [ 0.0, np.cos( eul[1] ), -np.sin( eul[1] ) ]
    C[2,:] = [ 0.0, np.sin( eul[1] ), np.cos( eul[1] ) ]

    B[0,:] = [ np.cos( eul[2] ), -np.sin( eul[2] ), 0.0 ] 
    B[1,:] = [ np.sin( eul[2] ), np.cos( eul[2] ), 0.0 ]
    B[2,:] = [ 0.0, 0.0, 1.0 ]

    R = np.matmul(D, C)
    R = np.matmul(R, B)

    return(R)

# -----------------------------------------------------------------------------

def bond(R):
    # From the euler rotation matrix, compute the 6-by-6 bond matrix
    M = np.zeros([6,6])
    M[0,:] = [ R[0,0]**2, R[0,1]**2, R[0,2]**2, 2*R[0,1]*R[0,2], 2*R[0,2]*R[0,0], 2*R[0,0]*R[0,1] ]
    M[1,:] = [ R[1,0]**2, R[1,1]**2, R[1,2]**2, 2*R[1,1]*R[1,2], 2*R[1,2]*R[1,0], 2*R[1,0] * R[1,1] ]
    M[2,:] = [ R[2,0]**2, R[2,1]**2, R[2,2]**2, 2*R[2,1]*R[2,2], 2*R[2,2]*R[2,0], 2*R[2,0] * R[2,1] ]
    M[3,:] = [ R[1,0]* R[2,0], R[1,1] * R[2,1], R[1,2] * R[2,2], R[1,1] * R[2,2] + R[1,2]*R[2,1], R[1,0]*R[2,2] + R[1,2]*R[2,0], R[1,1]*R[2,0] + R[1,0]*R[2,1] ]
    M[4,:] = [ R[2,0]* R[0,0], R[2,1] * R[0,1], R[2,2] * R[0,2], R[0,1] * R[2,2] + R[0,2]*R[2,1], R[0,2]*R[2,0] + R[0,0]*R[2,2], R[0,0]*R[2,1] + R[0,1]*R[2,0] ]
    M[5,:] = [ R[0,0]* R[1,0], R[0,1] * R[1,1], R[0,2] * R[1,2], R[0,1] * R[1,2] + R[0,2]*R[1,1], R[0,2]*R[1,0] + R[0,0]*R[1,2], R[0,0]*R[1,1] + R[0,1]*R[1,0] ]

    return(M)