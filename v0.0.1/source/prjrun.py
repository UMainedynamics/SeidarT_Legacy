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
    '-p', '--prjfile', nargs=1, type=str, required = True,
    help='the full file path for the project file', default=None
)

parser.add_argument(
    '-m', '--model', nargs = 1, type = str, required = False,
    help = """Specify whether to run the seismic (s), or electromagnetic (e), 
    or none (default = n)""",
    default = 'n'
)


# Get the arguments
args = parser.parse_args()
project_file = ''.join(args.prjfile)
model_type = ''.join(args.model)
pwd = os.path.dirname(project_file)

# ------------- Globals ----------------
clight = 2.99792458e8 # In general

# ============================ Create the objects =============================
# Let's initiate the domain
domain, material, seismic, electromag = loadproject(
    project_file,
    Domain(), 
    Material(),
    Model(),
    Model()
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

# Correct the domain to account for the absorbing boundary

# ---------------------------------------------------------------------
# We will always compute the coefficients but we need to make sure that we have
# everything needed to compute them
if seismic.exit_status == 0 and material.material_flag:
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

if electromag.exit_status == 0 and material.material_flag:
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


# ---------------------------------- SEISMIC ----------------------------------
# Do some seismic modeling
if model_type == 's':
    if seismic.exit_status == 0 and not seismic.compute_coefficients:
        # The coefficients are provided. We don't need the material
        
        print('Modeling the seismic wavefield.\n')
        seismic, domain = prepme(seismic, domain)
        
        # We need to set a density gradient at air interfaces because high
        # density gradients lead to numerical instability
        rhograd = airsurf(material, domain, 2)
        
        # Write the coefficient image to a fortran file
        seis25d.seismicfdtd25d.stiffness_write(
            domain.geometry+1,
            seismic.tensor_coefficients,
            domain.cpml,
            rhograd,
            domain.nx,
            domain.nz
        )
        if domain.dim == 2.5:
            print('Running 2.5D model')
            # Run the FDTD
            seis25d.seismicfdtd25d.seismic_cpml_25d(
                domain.nx + 2*domain.cpml, 
                domain.ny + 2*domain.cpml, 
                domain.nz + 2*domain.cpml,
                domain.dx,domain.dy, domain.dz,
                domain.cpml,
                seismic.src,
                seismic.f0,
                seismic.time_steps
            )
        else:
            seis2d.seismicfdtd2d.seismic_cpml_2d(
                domain.nx + 2*domain.cpml, 
                domain.nz + 2*domain.cpml,
                domain.dx, domain.dz,
                domain.cpml,
                seismic.src,
                seismic.f0,
                seismic.time_steps
            )

# ------------------------------ ELECTROMAGNETIC ------------------------------

# Now let's see if we can do some em modeling
if model_type == 'e':
    if electromag.exit_status == 0 and not electromag.compute_coefficients:
        # The coefficients are provided. We don't need the material values
        
        print('Modeling the electromagnetic wavefield.\n')
        # electromag, domain = prepme(electromag, domain)

        electromag, domain = prepme(electromag, domain)
        em25d.electromagfdtd25d.permittivity_write(
            domain.geometry+1,
            electromag.tensor_coefficients,
            domain.cpml,
            domain.nx, 
            domain.nz
        )
        
        if domain.dim == 2.5:
            print('Running 2.5D model')
            em25d.electromagfdtd25d.electromag_cpml_25d(
                domain.nx + 2*domain.cpml, 
                domain.ny + 2*domain.cpml, 
                domain.nz + 2*domain.cpml,
                domain.dx, domain.dy, domain.dz,
                domain.cpml,
                electromag.src,
                electromag.f0,
                electromag.time_steps,
            )
        else:
            print('Running 2 D model')
            em2d.electromagfdtd2d.electromag_cpml_2d(
                domain.nx + 2*domain.cpml,
                domain.nz + 2*domain.cpml,
                domain.dx, domain.dz,
                domain.cpml,
                electromag.src,
                electromag.f0,
                electromag.time_steps
            )
