#!/usr/bin/env python3

# This script will read the parameters given from a project file and run the 
# specified models. 

# -----------------------------------------------------------------------------

import argparse
import os.path
import numpy as np
import matplotlib.image as mpimg
from subprocess import call
import material_functions as mf
from class_definitions import *

# Modeling modules
import seismicfdtd2d as seis2d
import seismicfdtd25d as seis25d
import emfdtd2d as em2d
import emfdtd25d as em25d

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(
    description="""The SeidarT software requires a
    .PNG image that is used to construct the model domain for seismic and
    electromagnetic wave propagation. Given the image file, a project file
    will be constructed which contains all the necessary parameters to be
    read in to the finite differences time domain modeling schemes."""
)

parser.add_argument(
    'project_file', nargs=1, type=str,
    help='the full file path for the project file', default=None
)

parser.add_argument(
    '-M', '--model', nargs = 1, type = str, required = False,
    help = 'Specify whether to run the seismic (s), or electromagnetic (e), both (b) or none (default = n)',
    default = 'n'
)


# Get the arguments
args = parser.parse_args()
project_file = ''.join(args.project_file)
model_type = ''.join(args.model)
pwd = os.path.dirname(project_file)


# ============================ Create the objects =============================
# Let's initiate the domain
domain = Domain()
material = Material()
seismic = Model()
electromag = Model()


domain, material, seismic, electromag = loadproject(
    project_file,
    domain, material,
    seismic,
    electromag
)

# =================================== Model ===================================
# Check to make sure all the domain values were given
domain.para_check()

# Check the model inputs
seismic.para_check()
electromag.para_check()

# Check to make sure all the model inputs are satisfied
seismic.tensor_check()
electromag.tensor_check()

# Drop the rgb values 
material.material_list = np.delete(material.material_list, 2, axis = 1)
# print(material.material_list)

# Check the material list
material.para_check()

# ---------------------------------- SEISMIC ----------------------------------
# Let's see if we can do some seismic modeling
if model_type == 's' or model_type == 'b':
    if seismic.exit_status == 0 and not seismic.compute_coefficients:
        # The coefficients are provided. We don't need the material
        if model_type != 'n':
            seismic.time_steps = int(seismic.time_steps[0])
            seismic.f0 = float(seismic.f0[0])
            seismic.theta = float(seismic.theta[0])
            seismic.x = float(seismic.x[0])
            seismic.z = float(seismic.z[0])
            seismic.tensor_coefficients = seismic.tensor_coefficients.astype(float)
            
            print('Modeling the seismic wavefield.\n')
            
            if domain.dim == 2.5:
                # There are additional values we need to assign
                domain.ny = int(domain.ny[0])
                domain.dy = float(domain.dy[0])
                seismic.y = float(seismic.y[0])
                seismic.phi = float(seismic.phi[0])
                
                src = np.array(
                    [
                        seismic.x/domain.dx, 
                        seismic.y/domain.dy, 
                        seismic.z/domain.dz
                    ]
                ).astype(int)
                force=np.array([seismic.theta, seismic.phi])
                
                nx = domain.nx + 2*domain.cpml
                ny = domain.ny + 2*domain.cpml
                nz = domain.nz + 2*domain.cpml
                
                print('Running 2.5D model')
                print(seismic.time_steps)
                
                seis25d.seismicfdtd25d.stiffness_write(
                    domain.geometry+1,
                    seismic.tensor_coefficients,
                    domain.cpml,
                    domain.nx,
                    domain.nz
                )
                
                seis25d.seismicfdtd25d.seismic_cpml_25d(
                    nx, ny, nz,
                    domain.dx,
                    domain.dy,
                    domain.dz,
                    domain.cpml,
                    src,
                    seismic.f0,
                    seismic.time_steps
                )
            else:
                src = np.array([seismic.x/domain.dx, seismic.z/domain.dz]).astype(int)
                print('Running 2 D model')
                
                seis25d.seismicfdtd25d.stiffness_write(
                    domain.geometry+1,
                    seismic.tensor_coefficients,
                    domain.cpml,
                    domain.nx,
                    domain.nz
                )
                seis2d.seismicfdtd2d.seismic_cpml_2d(
                    nx, nz,
                    domain.dx,
                    domain.dz,
                    domain.cpml,
                    src,
                    seismic.f0,
                    seismic.time_steps,
                )

    elif seismic.exit_status == 0 and seismic.compute_coefficients and material.material_flag:
        # The coefficients aren't provided but the materials are so we can compute them
        # assign the materials to their respective corners
        
        print('Computing the stiffness coefficients.')
        material.sort_material_list()
        tensor = material.functions.get_seismic(
            temp = material.temp, rho = material.rho,
            pore = material.pore,wc = material.wc, abool = material.abool,
            angfile = material.angfiles, material_name = material.material
        )
        
        # Before we append the coefficients to the text file let's round to the second decimal
        tensor = np.round(tensor, 2)
        ind = np.where(tensor.max() == tensor)
        max_rho = tensor[ ind[0][0], -1]
        dt = np.min([domain.dx, domain.dz]) / np.sqrt(3.0 * tensor.max()/max_rho )
        
        # We're going to find the lines marked 'C' and input the values there
        append_coefficients(project_file, tensor, CP = 'C', dt = dt)
        print("Finished. Appending to project file.\n")
        
        if model_type != 'n':
            seismic.time_steps = int(seismic.time_steps[0])
            seismic.f0 = float(seismic.f0[0])
            seismic.theta = float(seismic.theta[0])
            seismic.x = float(seismic.x[0])
            seismic.z = float(seismic.z[0])
            seismic.tensor_coefficients = tensor
            
            print('Modeling the seismic wavefield.\n')
            
            if domain.dim == 2.5:
                # There are additional values we need to assign
                domain.ny = int(domain.ny[0])
                domain.dy = float(domain.dy[0])
                seismic.y = float(seismic.y[0])
                seismic.phi = float(seismic.phi[0])
                
                src = np.array(
                    [seismic.x/domain.dx, seismic.y/domain.dy, seismic.z/domain.dz]
                ).astype(int)
                force=np.array([seismic.theta, seismic.phi])
                
                nx = domain.nx + 2*domain.cpml
                ny = domain.ny + 2*domain.cpml
                nz = domain.nz + 2*domain.cpml 
                
                print('Running 2.5D model')
                seis25d.seismicfdtd25d.stiffness_write(
                    domain.geometry+1, 
                    seismic.tensor_coefficients,
                    domain.cpml,
                    domain.nx,
                    domain.nz
                )
                seis25d.seismicfdtd25d.seismic_cpml_25d(
                    nx, ny, nz,
                    domain.dx,
                    domain.dy,
                    domain.dz,
                    domain.cpml,
                    src,
                    seismic.f0,
                    seismic.time_steps
                )
            
            else:
                src = np.array([seismic.x/domain.dx, seismic.z/domain.dz]).astype(int)
                print('Running 2 D model')
                
                seis25d.seismicfdtd25d.stiffness_write(
                    domain.geometry+1,
                    seismic.tensor_coefficients,
                    domain.cpml,
                    domain.nx,
                    domain.nz
                )
                seis2d.seismicfdtd2d.seismic_cpml_2d(
                    nx, nz,
                    domain.dx,
                    domain.dz,
                    domain.cpml,
                    src,
                    seismic.f0,
                    seismic.time_steps,
                )
    else:
        # We don't have the materials or the coefficients
        print('Unable to model seismic. Check your project file for errors.\n')


# ------------------------------ ELECTROMAGNETIC ------------------------------
clight = 2.99792458e8 # In general

# Now let's see if we can do some em modeling
if model_type == 'e' or model_type == 'b':
    if electromag.exit_status == 0 and not electromag.compute_coefficients:
        # The coefficients are provided. We don't need the material
        if model_type != 'n':
            electromag.time_steps = int(electromag.time_steps[0])
            electromag.f0 = float(electromag.f0[0])
            electromag.theta = float(electromag.theta[0])
            electromag.x = float(electromag.x[0])
            electromag.z = float(electromag.z[0])
            electromag.tensor_coefficients = electromag.tensor_coefficients.astype(float)
            
            src = np.array([electromag.x/domain.dx, electromag.z/domain.dz]).astype(int)
            print('Modeling the electromagnetic wavefield.\n')
            
            if domain.dim == 2.5:
                # There are additional values we need to assign
                domain.ny = int(domain.ny[0])
                domain.dy = float(domain.dy[0])
                electromag.y = float(electromag.y[0])
                electromag.phi = float(electromag.phi[0])
                
                src = np.array(
                    [electromag.x/domain.dx, electromag.y/domain.dy, electromag.z/domain.dz]
                ).astype(int)
                force=np.array([electromag.theta, electromag.phi])
                nx = domain.nx + 2*domain.cpml
                ny = domain.ny + 2*domain.cpml
                nz = domain.nz + 2*domain.cpml 
                
                print('Running 2.5D model')
                
                em25d.electromagfdtd25d.permittivity_write(
                    domain.geometry+1,
                    electromag.tensor_coefficients,
                    domain.cpml,
                    domain.nx,
                    domain.nz
                )
                
                em25d.electromagfdtd25d.electromag_cpml_25d(
                    nx, ny, nz,
                    domain.dx,
                    domain.dy,
                    domain.dz,
                    domain.cpml,
                    src,
                    electromag.f0,
                    electromag.time_steps,
                    force
                )
            else:
                em2d.electromagfdtd2d.doall(
                    domain.geometry+1,
                    electromag.tensor_coefficients,
                    domain.dx,
                    domain.dz,
                    domain.cpml,
                    src,
                    electromag.f0,
                    electromag.time_steps,
                    electromag.theta
                )

    elif electromag.exit_status == 0 and electromag.compute_coefficients and material.material_flag:
        # The coefficients aren't provided but the materials are so we can compute them
        print('Computing the permittivity and conductivity coefficients.')
        material.sort_material_list()
        tensor = material.functions.get_perm(
            temp = material.temp,
            pore = material.pore,
            wc = material.wc,
            abool = material.abool,
            angfile = material.angfiles,
            material_name = material.material
        )
        dt = np.min([domain.dx, domain.dz])/(2.0*clight)
        
        append_coefficients(project_file, tensor, CP = 'P', dt = dt)
        print("Finished. Appending to project file.\n")
        
        if model_type != 'n':
            electromag.time_steps = int(electromag.time_steps[0])
            electromag.f0 = float(electromag.f0[0])
            electromag.theta = float(electromag.theta[0])
            electromag.x = float(electromag.x[0])
            electromag.z = float(electromag.z[0])
            electromag.tensor_coefficients = tensor
            
            src = np.array([electromag.x/domain.dx, electromag.z/domain.dz]).astype(int)
            print('Modeling the electromagnetic wavefield.\n')
            
            if domain.dim == 2.5:
                # There are additional values we need to assign
                domain.ny = int(domain.ny[0])
                domain.dy = float(domain.dy[0])
                electromag.y = float(electromag.y[0])
                electromag.phi = float(electromag.phi[0])
                
                src = np.array(
                    [electromag.x/domain.dx, electromag.y/domain.dy, electromag.z/domain.dz]
                ).astype(int)
                force=np.array([electromag.theta, electromag.phi])
                nx = domain.nx + 2*domain.cpml
                ny = domain.ny + 2*domain.cpml
                nz = domain.nz + 2*domain.cpml 
                
                print('Running 2.5D model')
                
                em25d.electromagfdtd25d.permittivity_write(
                    domain.geometry+1,
                    electromag.tensor_coefficients,
                    domain.cpml,
                    domain.nx,
                    domain.nz
                )
                em25d.electromagfdtd25d.electromag_cpml_25d(
                    nx, ny, nz,
                    domain.dx,
                    domain.dy,
                    domain.dz,
                    domain.cpml,
                    src,
                    electromag.f0,
                    electromag.time_steps,
                    force
                )
            else:
                em2d.electromagfdtd2d.doall(
                    domain.geometry+1,
                    electromag.tensor_coefficients,
                    domain.dx,
                    domain.dz,
                    domain.cpml,
                    src,
                    electromag.f0,
                    electromag.time_steps,
                    electromag.theta
                )
    else:
        # We don't have the materials neither the coefficients
        print('Unable to model electromagnetics. Check your project file for errors.')
