import numpy as np 
import material_functions as mf 
import matplotlib.image as mpimg
import os.path
from subprocess import call
from scipy.io import FortranFile
import glob

from imdefinitions import * 

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
    	self.write = None
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
    		self.dim = float(self.dim[0])
    		self.nx = int(self.nx[0])
    		self.nz = int(self.nz[0])
    		self.dx = float(self.dx[0])
    		self.dz = float(self.dz[0])
    		self.cpml = int(self.cpml[0])
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
		self.src = None
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


    # -------------------------- Function Definitions -------------------------
def getrcx(channel, rcx, domain):
    # input rcx as an n-by-2 array integer values for their indices. 
    all_files = glob.glob(channel + '*.dat')
    all_files.sort()
    m = len(all_files)
    n = len(rcx[:,1])
    timeseries = np.zeros([m,n])
    
    if domain.dim == '2':
        domain.ny = None
    
    if domain.dim == '2.5':
        for i, fn in enumerate(all_files, start = 0):
            npdat = read_dat(fn,channel,domain)
            for j in range(0, n):
                # Don't forget x is columns and z is rows
                timeseries[i,j] = npdat[
                    int(rcx[j,2]),
                    int(rcx[j,1]),
                    int(rcx[j,0])
				]
    else:
        for i, fn in enumerate(all_files, start = 0):
            npdat = read_dat(fn, channel, domain)
            for j in range(0, n):
                # Don't forget x is columns and z is rows
                timeseries[i,j] = npdat[
                    int(rcx[j,2]),
                    int(rcx[j,0])
                ]
    
    # Save the array as csv for other types of processing
    np.savetxt("receiver_array.csv", timeseries, delimiter = ",")


def read_dat(fn, channel, domain):
    if domain.dim == '2.5':
        if channel == 'Ex':
            NX = domain.nz
            NY = domain.ny
            NZ = domain.nx-1
        elif channel == 'Ey':
            NX = domain.nz
            NY = domain.ny-1
            NZ = domain.nx
        elif channel == 'Ez':
            NX = domain.nz-1
            NY = domain.ny
            NZ = domain.nx
        else:
            NX = domain.nz
            NY = domain.ny 
            NZ = domain.nx
    else:
        if channel == 'Ex':
        	NX = domain.nz
        	NZ = domain.nx-1
        elif channel == 'Ez':
            NX = domain.nz-1
            NZ = domain.nx
        else:
            NX = domain.nz 
            NZ = domain.nx
    
    
    f = FortranFile(fn, 'r')
    dat = f.read_reals(dtype = 'float32')
    
    if domain.dim == '2.5':
        dat = dat.reshape(NX, NY, NZ)
    else:
        dat = dat.reshape(NX, NZ)
    
    f.close()
    return(dat)

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
# ========================= Read/Assign Project File ==========================
# =============================================================================
def loadproject(project_file, domain, material, seismic, electromag):
    # domain, material, seismic, electromag are the objects that we will assign
    # values to
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
			if temp[1] == 'write':
				domain.write = temp[2].rsplit()


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
				seismic.tensor_coefficients = np.row_stack((seismic.tensor_coefficients, temp[1:]))

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
	return domain, material, seismic, electromag


# -----------------------------------------------------------------------------
# Make sure variables are in the correct type for Fortran
def prepme(modobj, domain):
    # Check if there are no errors and the coefficients have been computed
    modobj.time_steps = int(modobj.time_steps[0])
    modobj.f0 = float(modobj.f0[0])
    modobj.theta = float(modobj.theta[0])
    modobj.x = float(modobj.x[0])
    modobj.z = float(modobj.z[0])
    modobj.tensor_coefficients = modobj.tensor_coefficients.astype(float)
    # Put source and domain parameters in correct order
    if domain.dim == 2.5:
        # There are additional values we need to assign
        domain.ny = int(domain.ny[0])
        domain.dy = float(domain.dy[0])
        modobj.y = float(modobj.y[0])
        modobj.phi = float(modobj.phi[0])
        
        modobj.src = np.array(
            [
            	modobj.x/domain.dx, 
            	modobj.y/domain.dy, 
            	modobj.z/domain.dz
        	]
        ).astype(int)
    else:
        modobj.src = np.array(
            [
                modobj.x/domain.dx, 
                modobj.z/domain.dz
            ]
        ).astype(int)
    
    return(modobj, domain)

# ----------------------
# Append coefficients
def coefs2prj(modobj, matobj, domobj, modtype):
    pass


def airsurf(material, domain, N = 2):
    # This can be generalized a little better, but for now...
    airnum = material.material_list[material.material_list[:,1] == 'air', 0]
    
    if airnum:
        airnum = int(airnum[0])
        gradmatrix = (domain.geometry != airnum).astype(int)
        # Take the gradient in both directions
        gradz = np.diff(gradmatrix, axis = 0)
        gradx = np.diff(gradmatrix, axis = 1)
        
        # For grady we will append a column of zeros at the beginning so that the value
        # 1 is located at the interface but on the non-air side 
        gradzpos = np.row_stack([np.zeros([gradz.shape[1] ]),gradz])
        # For gradx we will append a row of zeros at the beginning 
        gradxpos = np.column_stack([np.zeros([gradx.shape[0] ]),gradx]) 
        # -1 also means that there is an air interface. We will need to append to the
        # end of the array, then we can just flip the sign
        gradzneg = (np.row_stack( [gradz, np.zeros([gradz.shape[1]]) ] ) )
        gradxneg = (np.column_stack( [gradx, np.zeros([gradx.shape[0]]) ] ) )
        
        # At the surface we want to have 15% of density
        grad = np.zeros( [gradx.shape[0], gradz.shape[1], N] ) 
        grad[:,:,0] = gradzpos + gradxpos - gradzneg - gradxneg
        grad[ grad[:,:,0]>0, 0] = 0.15 
        
        # We will make the change gradational by splitting the difference each step 
        # For instance, 1 - 0.15 = 0.85, so the next percentage will be 
        # 0.85/2 + 0.15 and so on
        pct = np.zeros([N])
        pct[0] = 0.15
        
        for ind in range(1, N):
            pct[ind] = pct[ind-1] + (1-pct[ind-1])/2
            gradzpos = np.row_stack( [np.zeros([gradz.shape[1] ]),gradzpos] )[:-1,:]
            gradxpos = np.column_stack( [np.zeros( [ gradx.shape[0] ] ),gradxpos ])[:,:-1]
            
            gradzneg = (np.row_stack( [ gradzneg, np.zeros( [gradz.shape[1] ]) ] )[1:,:])
            gradxneg = (np.column_stack( [ gradxneg, np.zeros( [gradx.shape[0] ]) ] )[:,1:])
            grad[:,:, ind] = gradzpos + gradxpos - gradzneg - gradxneg
            grad[ grad[:,:,ind] > 0, ind] = pct[ind]
        
        gradcomp = np.zeros( [grad.shape[0], grad.shape[1] ])
        for ind in range(N-1, -1, -1):
            gradcomp[ grad[:,:,ind] > 0] = pct[ind]
        
        gradcomp[ gradcomp == 0] = 1
    else:
        gradcomp = np.ones([int(domain.nx), int(domain.nz) ])
    
    return(gradcomp)

def indvar(N, modobj, domain):
    x = np.linspace(
        0, 
        float(domain.dx[0]) * (nx - 1), 
        N[0]
    )
    z = np.linspace(
        0, 
        float(domain.dz[0]) * (nz - 1), 
        N[2]
    )
    t = np.linspace(
		0,
  		float( modobj.dt[0]) * (float(modobj.time_steps[0]) - 1),
		int( modobj.time_steps[0] )
	)
    try:
		y = np.linspace(
    		0, 
        	float(domain.dy[0]) * (ny - 1), 
        	N[1]
    	)
	except:
		y = None
	
    return(x,y,z,t)
    