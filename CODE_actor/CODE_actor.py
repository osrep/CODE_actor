# -*- coding: utf-8 -*-
"""
Created on Tue May 3 12:00:00 2020

@author: Soma Olasz
"""

# Import necessary modules and Python scripts.
import hdf5Write

import numpy as np
import matlab.engine
import itertools
import pickle
import ual
from time import asctime
import copy


# Start a Matlab engine to be able to call Matlab scripts
eng = matlab.engine.start_matlab()

try:
	
	# Read in data from CPO
	fBig = distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].values.scalar[:]
	exty = distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[0].objects[0].geo[:,0,0,0]
	ext_rho_tor = distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[1].objects[0].geo[:,0,0,0]

	externalSwitch = 1
	
except IndexError:

	externalSwitch = 0

if externalSwitch == 1:

	ext_rho_size = ext_rho_tor.size
	
	# Extract numerical parameters necessary to recreate the grid of the iput distribution

	yMax = max(exty)
	Ny = exty.size
	
	Nxi = fBig.size/ext_rho_size/Ny

	# Calulate the matrix size, used for the reshaping of the external distribution
	matrixSize = Nxi*Ny

	eng.workspace["yMax"] = float(yMax)
	eng.workspace["Ny"] = float(Ny)
	eng.workspace["Nxi"] = float(Nxi)

	# Extract the distribution for the different rho coordinates
	f = np.zeros((ext_rho_size,matrixSize))
	
	for i in range (0,ext_rho_size):

		f[i,:] = fBig[0:matrixSize]

		fBig = fBig[matrixSize:]

	
	for i in range (0, ext_rho_size):
		
		# take the part of the distribution function for the current rho coordinate.
		temp = f[i, :]
		temp = temp.tolist()
		
		# If first rho indey, save it as the distribution
		if i == 0:
			eng.workspace["fTotal"] = temp
			eng.eval("fTotal = cell2mat(transpose(fTotal));", nargout = 0)
		
		# Else, append it to the existing distribution array
		else:
			eng.workspace["temp"] = temp
			eng.eval("temp = cell2mat(transpose(temp));", nargout = 0)
			eng.eval("fTotal = cat(2, fTotal, temp);", nargout = 0)


else:
	Ny = 60
	Nxi = 10
	yMax = 200
	# get the grid parameters from the input grid
	eng.workspace["yMax"] = float(yMax)
	eng.workspace["Ny"] = float(Ny)
	eng.workspace["Nxi"] = float(Nxi)


# Set the working directory to the CODE directory
eng.cd('/pfs/work/g2solasz/git/CODE_actor/ext/CODE/src')
eng.addpath('/pfs/work/g2solasz/git/CODE_actor/CODE_actor')

eng.InitCODE(nargout=0)

# Create the plasma state object. Give normalization reference values for T and n.
TRef = 10.0;
nRef = 1e19;

eng.workspace["TRef"] = float(TRef)
eng.workspace["nRef"] = float(nRef)

eng.eval("S = State(TRef,nRef);", nargout = 0)

## Set the momentum parameters
eng.workspace["yGridMode"] = float(4)
eng.eval("S.momentumGrid.setyMax(yMax);", nargout = 0)
eng.eval("S.momentumGrid.setNy(Ny);", nargout = 0)
eng.eval("S.momentumGrid.setNxi(Nxi);", nargout = 0)
eng.eval("S.momentumGrid.setyGridMode(yGridMode);", nargout = 0)


## Set the TimeGrid
# Read in the workflow time step size
tMax = parameters["dt_in"] # in s
dt = 1
timeStepMode = 0

eng.workspace["tMax_s"] = float(tMax)
eng.eval("tMax = tMax_s*S.reference.nueeRef;", nargout = 0) # convert to normalized units
eng.workspace["dt"] = float(dt)
eng.workspace["timeStepMode"] = float(timeStepMode)
eng.eval("S.timeGrid.settMax(tMax);", nargout = 0)
eng.eval("S.timeGrid.setdt(dt);", nargout = 0)
eng.eval("S.timeGrid.settimeStepMode(timeStepMode);", nargout = 0)

## Set Physical Quantities from CPOs
# Get the number of rho coordinates
rho_size = size(coreprof0[0].rho_tor_norm)

# Initialize arrays for physical parameters
temperature = np.zeros(rho_size)
density = np.zeros(rho_size)
Z_eff = np.zeros(rho_size)
B0 = np.zeros(rho_size)
E_parallel = np.zeros(rho_size)


rhoTor_arr = np.zeros(rho_size)
time = coreprof0[0].time

# Initialize variables for storing calculation results
totalDistribution = []		# data storage for CPOs
yFinal = []			# data storage for CPOs

# Constant physical parameters
for i in range(rho_size):
	
	# Fill physics arrays with values from CPOs
	temperature[i] = coreprof0[0].te.value[i]					# eV
	density[i] = coreprof0[0].ne.value[i]						# m^{-3}
	Z_eff[i] = coreprof0[0].profiles1d.zeff.value[i]				# Z_eff
	B0[i] = coreprof0[0].toroid_field.b0						# T
	E_parallel[i] = coreprof0[0].profiles1d.eparallel.value[i]			# V/m

	rhoTor_arr[i] = coreprof0[0].rho_tor[i]						# m

# Prepare the CODE for the various rho values and run the calculation
for i in range(rho_size):
	
	eng.workspace["T"] = float(temperature[i])
	eng.workspace["n"] = float(density[i])
	eng.workspace["Z"] = float(Z_eff[i])
	eng.workspace["E"] = float(E_parallel[i])
	eng.workspace["B"] = float(B0[i])

	eng.eval("S.physicalParams.setPhysicalParameters(T,n,Z,E,B);", nargout = 0)

	## Set the Operator settings for the run
	# Create Equation settings based on the set plasma parameters
	eng.eval("es = EquationSettings(S);", nargout = 0)

	# Various operator settings can be made, see EquationSettings.m for details 
	eng.eval("es.setsourceMode(0);", nargout = 0)
	eng.eval("es.setcollisionOperator(0);", nargout = 0)
	eng.eval("es.setbremsMode(0);", nargout = 0)

	## Create the Solver object based on State and equation Settings
	eng.eval("solver = SmartLUSolver(S,es);", nargout = 0)
	
	if externalSwitch == 1:
		
		# Extract the distribution function for the current rho index
		eng.workspace["i"] = i + 1
		eng.eval("f = fTotal(:, i);", nargout = 0)
		
		# Create a momentum grid object for the input distribution based on the State reference
		eng.eval("inputMomentumGrid = MomentumGrid(S.reference, 'Ny', Ny, 'yMax', yMax, 'Nxi', Nxi);", nargout = 0)

		# Create a distribution object for the input distribution using the State object
		eng.eval("inputDistribution = Distribution(f, S, 1);", nargout = 0)
		
		# Set the distribution
		eng.eval("solver.setf(inputDistribution, inputMomentumGrid);", nargout = 0)
	
	else:
		eng.eval("solver.setftoMaxwellian();", nargout = 0)

	# Set the number of time steps to be saved
	eng.eval("solver.setnSaves(50);", nargout = 0)

	# Perform calculation
	eng.eval("output = solver.takeTimeSteps();", nargout = 0)
	
	output_f = eng.eval("output.distributions{end}.f;")
	output_y = eng.eval("output.distributions{end}.momentumGrid.y;")
	
	# Take the data from the output, which will go into the CPO.
	output_f = np.array(output_f).tolist()
	output_y = np.array(output_y).tolist()
	
	# flatten the data to python list, so it can be given to the CPO
	output_f = list(itertools.chain.from_iterable(output_f))
	output_y = list(itertools.chain.from_iterable(output_y))
	
	
	# Save coordinates from first calculation to write to CPO
	if i == 0:
		yFinal = output_y
	
	# Check if the grids are the same for the later calculations as the saved grid
	else:
		if not  yFinal == output_y:
			raise Exception('The y grid is not the same for rho index {} as for the first index.'.format(i+1))
		
	totalDistribution += output_f

# Convert rho coordinates to list so it can be given to CPOs
rhoTor = rhoTor_arr.tolist()

# Write calculation results to CPO
# Give run and shot numbers needed for CPO writing
shot = parameters["shotnumber"]
run = parameters["run_out"]

# Initialize CPO structure
distribution0 = ual.distribution.distributionArray()

# Set numerical parameters for CPO writing
nCoord = 2	# number of different coordinates used (y, rho_tor)- Not sure what to do with Xi TODO

# initialize the CPO
distribution0.resize(1)
distribution0.array[0].distri_vec.resize(1)
distribution0.array[0].distri_vec[0].dist_func.f_expansion.resize(1)
distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces.resize(nCoord)

# fill the coordinates
for i in range (nCoord):
	distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects.resize(1)
	distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].coordtype.resize(1,1)

	# y coordinate
	if i == 0: 
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects[0].geo.resize(Ny,1,1,1)
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].coordtype[0,0] = 123
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects[0].geo[:,0,0,0] = yFinal
		
	# Xi cordinate should be given here, but avoided to prevent misinterpretation. It is not given as it is Legendre polinomial coordinates, which cannot be given to CPOS.
	
	# rho coordinate
	else:
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects[0].geo.resize(rho_size,1,1,1)
		# 107 is the coordinate convention for rho_tor (see https://portal.eufus.eu/documentation/ITM/html/itm_enum_types__coordinate_identifier.html#itm_enum_types__coordinate_identifier). Might have to change this later.
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].coordtype[0,0] = 107
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects[0].geo[:,0,0,0] = rhoTor
			
# Write the distribution to the CPO
distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].values.scalar.resize(Ny*Nxi*rho_size)
distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].values.scalar[:] = totalDistribution
		
# Write the time
distribution0.array[0].time = time


# Reshape data for hdf5 writing
temperature = temperature.reshape(1,rho_size)
density = density.reshape(1,rho_size)
#EHat = EHat.reshape(1,rho_size)				TODO not saved yet
Z_eff = Z_eff.reshape(1,rho_size)
B0 = B0.reshape(1,rho_size)
rhoTor = rhoTor_arr.reshape(1,rho_size)
E_parallel = E_parallel.reshape(1,rho_size)
#E_critical = E_critical.reshape(1,rho_size)			TODO not saved yet
time = np.array(time).reshape(1,1)

Distribution = np.array([np.transpose(Distribution[0,:,:])])


# Create dictionaries of the hdf5 input parameters
temperature = {"Name": 'temperature', "Data": temperature}
density = {"Name": 'density', "Data":  density}
#EHat = {"Name": 'EHat', "Data":  EHat}				TODO not saved yet
Z_eff = {"Name": 'Z_eff', "Data":  Z_eff}
B0 = {"Name": 'B0', "Data":  B0}
rhoTor = {"Name": 'rhoTor', "Data":  rhoTor}
E_parallel = {"Name": 'E_parallel', "Data":  E_parallel}
#E_critical = {"Name": 'E_critical', "Data":  E_critical}	TODO not saved yet
time = {"Name": 'time', "Data":  time}
Distribution = {"Name": 'Distribution', "Data":  Distribution}

# Put dictionaries into a list
hdf5_param_data = [temperature, density, EHat, Z_eff, B0, rhoTor, E_parallel, E_critical, time]
hdf5_dist_data = [Distribution]

# Write data to hdf5 file
hdf5Write.write_params(shot, run, hdf5_param_data)
hdf5Write.write_dist(shot, run, hdf5_dist_data)
