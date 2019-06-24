#!/usr/bin/env python3
#
# This script will read the parameters given from a project file and run the 
# specified models. 
#
# -----------------------------------------------------------------------------

import argparse
import os.path
import numpy as np
import matplotlib.image as mpimg
from subprocess import call
import material_functions as mf

# Modeling modules
import seismicfdtd2d as seis2d
import emfdtd2d as em2d

# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(description="""The seisarT software requires a 
	.PNG image that is used to construct the model domain for seismic and 
	electromagnetic wave propagation. Given the image file, a project file 
	will be constructed which contains all the necessary parameters to be 
	read in to the finite differences time domain modeling schemes.""" )

parser.add_argument( 'project_file', nargs=1, type=str, 
						help='the full file path for the project file', default=None)

parser.add_argument( '-M', '--model', nargs = 1, type = str, required = False,
	help = 'Specify whether to run the seismic (s), or electromagnetic (e), both (b) or none (default = n)',
	default = 'n')



# Get the arguments
args = parser.parse_args()
project_file = ''.join(args.project_file)
model_type = ''.join(args.model)
pwd = os.path.dirname(project_file)


# =============================================================================
# =========================== Define Class Variables ==========================
# =============================================================================

# We need to define some class variables
class Domain:
    def __init__(self):
        super().__init__()
        self.build()

    def build(self):
    	# Initialize variables
    	self.geometry = None
    	self.dim = None
    	self.nx = None
    	self.ny = None
    	self.nz = None
    	self.dx = None
    	self.dy = None
    	self.dz = None
    	self.cpml = None
    	self.exit_status = 1

    	# Flag to specify whether the all inputs are fulfilled
    	self.seismic_model = False 
    	self.electromag_model = False

    def para_check(self):
    	self.exit_status = 0
    	# make sure there's a geometry. This implies whether nx, and nz exist 
    	if self.geometry is None:
    		self.exit_status = 1
    		print('No geometry loaded\n')
    	if not self.dx or self.dx is None:
    		self.exit_status = 1
    		print('No step size in the x-direction')
    	if not self.dz or self.dz is None:
    		self.exit_status = 1
    		print('No step size in the y-direction')
    	if not self.cpml or self.cpml is None:
    		self.exit_status = 1
    		print('No cpml thickness is given')

    	# Check 2.5D
    	if self.dim == '2.5' and exit_status == 0: 
    		if self.ny is None or self.ny == 'n/a':
    			print('No dimension given for y-direction. Assigning default ny=3.')
    			self.ny = 3

    		if self.dy is None or self.dy == 'n/a':
    			self.dy = np.min( [int(self.dx), int(self.dz)])
    			print('No step size given for y-direction. Assigning min(dx,dz).')
    	# Convert variables to required types. Ny and Dy are already converted if given
    	if self.exit_status == 0:
    		domain.dim = float(domain.dim[0])
    		domain.nx = int(domain.nx[0])
    		domain.nz = int(domain.nz[0])
    		domain.dx = float(domain.dx[0])
    		domain.dz = float(domain.dz[0])
    		domain.cpml = int(domain.cpml[0])
    	else:
    		print('\n Domain inputs are not satisfied. I can"t go on anymore. \n')
    		# quit()
	
# -----------------------------------------------------------------------------

class Material:
    # initialize the class
    def __init__(self):
        super().__init__()
        self.build()

    def build(self):
    	self.material_list = np.array([]) # initialize
    	self.material_flag = False # Whether the materials were read in
    	
    	# We will assign each of the list variables
    	self.material = None
    	self.temp = None
    	self.attenuation = None
    	self.rho = None
    	self.pore = None
    	self.wc = None    	
    	self.abool = None
    	self.angfiles = None

    	# The processing functions
    	self.functions = mf

    def sort_material_list(self):
    	self.material = self.material_list[:,1]
    	self.temp = self.material_list[:,2].astype(float)
    	self.attenuation = self.material_list[:,3].astype(float)
    	self.rho = self.material_list[:,4].astype(float)
    	self.pore = self.material_list[:,5].astype(float)
    	self.wc = self.material_list[:,6].astype(float)
    	self.abool = self.material_list[:,7] == 'True'
    	self.angfiles = self.material_list[:,8]

    def para_check(self):
    	# The fields for the materials in the input are defined as: 
    	# 'id, R/G/B, Temp. Atten., Dens., Por., WC, Anis, ANG_File'
    	# but the R/G/B column is deleted

    	if len(self.material_list) > 0:
    		# Check to make sure the necessary fields are provided
    		check = 0

    		for row in self.material_list[:,0:7]:
    			for val in row:
    				if not val:
    					check = check + 1

    	if check == 0:
	    	file_check = 0
	    	for ind in range(0, self.material_list.shape[0]):
	    		if self.material_list[ind,7] == 'True' and not self.material_list[ind,8] or self.material_list[ind,8] == 'n/a':
	    			file_check = file_check + 1
    	else:
    		print('Material inputs aren"t satisfied.')
    		# quit()

    	if check == 0:
    		if file_check == 0:
    			self.material_flag = True
    		else:
    			print('No .ANG file specified for anisotropic material')

# -----------------------------------------------------------------------------
class Model:
	def __init__(self):
		super().__init__()
		self.build()

	def build(self):
		self.dt = None
		self.time_steps = None
		self.x = None
		self.y = None 
		self.z = None
		self.f0 = None
		self.theta = None
		self.phi = None
		self.tensor_coefficients = None
		self.compute_coefficients = True
		self.exit_status = 0

	def tensor_check(self):

		# If the tensors are there 
		check = 0

		# if not self.tensor_coefficients:
		# 	print('ldkjf')

		for row in self.tensor_coefficients:
			for val in row:
				if not val:
					check = check + 1

		if check == 0:
			self.compute_coefficients = False

	def para_check(self):

		if not self.time_steps:
			self.exit_status = 1
			print('Number of time steps aren"t satisfied.')

		if not self.x or not self.z:
			self.exit_status = 1
			print('No source location is specified.')

		if not self.f0:
			self.exit_status = 1
			print('No source frequency is specified.')

		# in case theta or phi aren't specified we can assign defaults
		if not self.theta:
			self.theta = 0

		# if y is specified but not phi
		if self.y and not self.phi:
			self.phi = 0


# -----------------------------------------------------------------------------



# =============================================================================
# ============================== Useful Functions =============================
# =============================================================================

def image2int(imfilename):
    # read the image
    img = mpimg.imread(imfilename)

    # Convert RGB to a single value
    rgb_int = np.array(65536*img[:,:,0] +  255*img[:,:,1] + img[:,:,2])
    
    # Get the unique values of the image
    rgb_uni = np.unique(rgb_int)
    
    # mat_id = np.array( range(0, len(rgb_uni) ) )
    for ind in range(0, len(rgb_uni) ):
    	rgb_int[ rgb_int == rgb_uni[ind] ] = ind	

    return rgb_int.astype(int)

# -----------------------------------------------------------------------------
# After computing tensor coefficients we want to append them to the given text
# file. In order to do this, we need to create a new file then move it to the 
# original file

def append_coefficients(prjfile, tensor, CP = None, dt=1 ):
	# CP is the line identifier (C - stiffness, P - permittivity). This has the
	# ability to be optional since there will is a difference between a 2nd 
	# order tensor and a 4th order tensor in regards to length but we might
	# want to include other types of modeling in the future.
	newfile = 'newfile.txt'

	ogprj = open(prjfile, 'r')
	temp = open(newfile, 'a')
	# We need to append the dt for each modeling type
	if CP == 'C':
		mt = 'S'
	else:
		mt = 'E'

	for line in ogprj.readlines():
		if line[0] == CP:
			line = line.split(',')
			temp.write( CP + ',' + ','.join( tensor[ int(float(line[1])),:].astype(str) ) + '\n' )
		elif line[0] == mt and line[2:4] == 'dt':
			temp.write( mt + ',dt,' + str(dt) + '\n' )
		else:
			temp.write(line)

		# if line[0] == mt:
		# 	line = line.split(',')
		# 	if line[1] == 'dt':
		# 		temp.write( mt + ',dt,' + str(dt) + '\n' )



	call('mv ' + newfile + ' ' + prjfile, shell = True)



# =============================================================================
# ============================ Create the objects =============================
# =============================================================================

# Let's initiate the domain
domain = Domain()
material = Material()
seismic = Model()
electromag = Model()



# =============================================================================
# ========================= Read/Assign Project File ==========================
# =============================================================================

f = open(project_file)

# Let the if train begin
for line in f:

	if line[0] == 'I':
		# There's a trailing new line value
		im = image2int(line[2:-1])
		domain.geometry = im.transpose().astype(int)

	# All domain inputs must be input except for ny and dy
	if line[0] == 'D':
		temp = line.split(',')

		if temp[1] == 'dim':
			domain.dim = temp[2].rsplit()
		if temp[1] == 'nx':
			domain.nx = temp[2].rsplit()
		if temp[1] == 'ny':
			domain.ny = temp[2].rsplit()
		if temp[1] == 'nz':
			domain.nz = temp[2].rsplit()
		if temp[1] == 'dx':
			domain.dx = temp[2].rsplit()
		if temp[1] == 'dy':
			domain.dy = temp[2].rsplit()
		if temp[1] == 'dz':
			domain.dz = temp[2].rsplit()
		if temp[1] == 'cpml':
			domain.cpml = temp[2].rsplit()

	if line[0] == 'S':
		temp = line.split(',')
		if temp[1] == 'dt':
			seismic.dt = temp[2].rsplit()
		if temp[1] == 'time_steps':
			seismic.time_steps = temp[2].rsplit()
		if temp[1] == 'x':
			seismic.x = temp[2].rsplit()
		if temp[1] == 'y':
			seismic.y = temp[2].rsplit()
		if temp[1] == 'z':
			seismic.z = temp[2].rsplit()
		if temp[1] == 'f0':
			seismic.f0 = temp[2].rsplit()
		if temp[1] == 'theta':
			seismic.theta = temp[2].rsplit()
		if temp[1] == 'phi':
			seismic.phi = temp[2].rsplit()

	if line[0] == 'E':
		temp = line.split(',')
		if temp[1] == 'dt':
			electromag.dt = temp[2].rsplit()
		if temp[1] == 'time_steps':
			electromag.time_steps = temp[2].rsplit()
		if temp[1] == 'x':
			electromag.x = temp[2].rsplit()
		if temp[1] == 'y':
			electromag.y = temp[2].rsplit()
		if temp[1] == 'z':
			electromag.z = temp[2].rsplit()
		if temp[1] == 'f0':
			electromag.f0 = temp[2].rsplit()
		if temp[1] == 'theta':
			electromag.theta = temp[2].rsplit()
		if temp[1] == 'phi':
			electromag.phi = temp[2].rsplit()
		
	if line[0] == 'M':
		line = line[0:-1]
		temp = line.split(',')
		if temp[1] == '0':
			material.material_list = temp[1:]
		else:
			material.material_list = np.row_stack( (material.material_list, temp[1:]))

	if line[0] == 'C':
		temp = line.split(',')
		# We need to remove the '\n' at the end. Whether the coefficients are 
		# given results in a different string
		try:
			temp[-1] = temp[-1].rsplit()[0]
		except:
			temp[-1] = ''

		if temp[1] == '0' or temp[1] == '0.0':
			seismic.tensor_coefficients = temp[1:]
		else:
			seismic.tensor_coefficients = np.row_stack( (seismic.tensor_coefficients, temp[1:]))

	if line[0] == 'P':
		temp = line.split(',')
		try:
			temp[-1] = temp[-1].rsplit()[0] # An index error will be given if coefficients are provided
		except:
			temp[-1] = ''

		if temp[1] == '0' or temp[1] == '0.0':
			electromag.tensor_coefficients = temp[1:]
		else:
			electromag.tensor_coefficients = np.row_stack( (electromag.tensor_coefficients, temp[1:]))

f.close()


# =============================================================================
# =================================== Model ===================================
# =============================================================================

# Check to make sure all the domain values were given
domain.para_check()

# check the model inputs
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

# Let's see if we can do some seismic modeling
if seismic.exit_status == 0 and not seismic.compute_coefficients:
	# The coefficients are provided. We don't need the material

	if model_type == 's' or model_type == 'b':
		seismic.time_steps = int(seismic.time_steps[0])
		seismic.f0 = float(seismic.f0[0])
		seismic.theta = float(seismic.theta[0])
		seismic.x = float(seismic.x[0])
		seismic.z = float(seismic.z[0])
		seismic.tensor_coefficients = seismic.tensor_coefficients.astype(float)

		src = np.array([seismic.x/domain.dx, seismic.z/domain.dz]).astype(int)
		print('Modeling the seismic wavefield.\n')

		seis2d.seismicfdtd2d.doall(domain.geometry+1, seismic.tensor_coefficients, 
			domain.dx, domain.dz, domain.cpml, src, seismic.f0, 
			seismic.time_steps, seismic.theta)

elif seismic.exit_status == 0 and seismic.compute_coefficients and material.material_flag:
	# The coefficients aren't provided but the materials are so we can compute them
	# assign the materials to their respective corners
	
	print('Computing the stiffness coefficients.')
	material.sort_material_list()
	tensor = material.functions.get_seismic(temp = material.temp, rho = material.rho, 
		pore = material.pore,wc = material.wc, abool = material.abool, 
		angfile = material.angfiles, material_name = material.material)

	# Before we append the coefficients to the text file let's round to the second decimal
	tensor = np.round(tensor, 2)

	ind = np.where(tensor.max() == tensor)
	max_rho = tensor[ ind[0][0], -1]
	dt = np.min([domain.dx, domain.dz])/(2.0*np.sqrt(tensor.max()/max_rho ))

	# We're going to find the lines marked 'C' and input the values there
	append_coefficients(project_file, tensor, CP = 'C', dt = dt)
	print("Finished. Appending to project file.\n")
	
	if model_type == 's' or model_type == 'b':
		seismic.time_steps = int(seismic.time_steps[0])
		seismic.f0 = float(seismic.f0[0])
		seismic.theta = float(seismic.theta[0])
		seismic.x = float(seismic.x[0])
		seismic.z = float(seismic.z[0])
		seismic.tensor_coefficients = tensor

		src = np.array([seismic.x/domain.dx, seismic.z/domain.dz]).astype(int)
		print('Modeling the seismic wavefield.\n')

		seis2d.seismicfdtd2d.doall(domain.geometry+1, seismic.tensor_coefficients, 
			domain.dx, domain.dz, domain.cpml, src, seismic.f0, 
			seismic.time_steps, seismic.theta)
else:
	# We don't have the materials or the coefficients
	print('Unable to model seismic. Check your project file for errors.\n')

# -----------------------------------------------------------------------------

clight = 2.99792458e8


# Now let's see if we can do some em modeling	
if electromag.exit_status == 0 and not electromag.compute_coefficients:
	# The coefficients are provided. We don't need the material

	if model_type == 'e' or model_type == 'b':

		electromag.time_steps = int(electromag.time_steps[0])
		electromag.f0 = float(electromag.f0[0])
		electromag.theta = float(electromag.theta[0])
		electromag.x = float(electromag.x[0])
		electromag.z = float(electromag.z[0])
		electromag.tensor_coefficients = electromag.tensor_coefficients.astype(float)

		src = np.array([electromag.x/domain.dx, electromag.z/domain.dz]).astype(int)
		print('Modeling the electromagnetic wavefield.\n')

		em2d.electromagfdtd2d.doall(domain.geometry+1, electromag.tensor_coefficients, 
			domain.dx, domain.dz, domain.cpml, src, electromag.f0, 
			electromag.time_steps, electromag.theta)

elif electromag.exit_status == 0 and electromag.compute_coefficients and material.material_flag:
	
	# The coefficients aren't provided but the materials are so we can compute them
	print('Computing the permittivity and conductivity coefficients.')
	material.sort_material_list()
	tensor = material.functions.get_perm(temp = material.temp, pore = material.pore,
		wc = material.wc, abool = material.abool, angfile = material.angfiles, 
		material_name = material.material)
	dt = np.min([domain.dx, domain.dz])/(2.0*clight)

	append_coefficients(project_file, tensor, CP = 'P', dt = dt)
	print("Finished. Appending to project file.\n")

	if model_type == 'e' or model_type == 'b':

		electromag.time_steps = int(electromag.time_steps[0])
		electromag.f0 = float(electromag.f0[0])
		electromag.theta = float(electromag.theta[0])
		electromag.x = float(electromag.x[0])
		electromag.z = float(electromag.z[0])
		electromag.tensor_coefficients = tensor

		src = np.array([electromag.x/domain.dx, electromag.z/domain.dz]).astype(int)
		print('Modeling the electromagnetic wavefield.\n')
		
		em2d.electromagfdtd2d.doall(domain.geometry+1, electromag.tensor_coefficients, 
			domain.dx, domain.dz, domain.cpml, src, electromag.f0, 
			electromag.time_steps, electromag.theta)

else:
	# We don't have the materials neither the coefficients
	print('Unable to model electromagnetics. Check your project file for errors.')
