# ########################### Main Home Simulator #################################################
#
# This module contains the main code of the Low Voltage Network Simulator (LVNS). The LVNS is an 
# application developed in python 2.7, based on the LVNS platform proposed by Torquato et al. 
# in: "A Monte Carlo Simulation Platform for Studying Low Voltage Residential Networks". The LVNS 
# can be used for modelling and simulating low voltage residential networks. The LVNS generates a 
# quasi-steady-state/harmonic network condition over an extended period such as 24 h. It consists 
# of two main components. A multiphase network model for power flow, harmonic, and motor starting 
# study capabilities, and the characteristics of various loads and generators based on time-of-use 
# probability curves.
# In order to launch the LVNS platform, the "main_home_simulator" model has to be called using 
# the following statement:
#	>> python main_home_simulator.py <input_file.txt>
# Where, input_file.txt is the input data file that contains the houses and electrical appliances 
# information (See the example in 3houses_base.txt)
#
# #################################################################################################
import sys
# #################################################################################################
# Input arguments validation
# #################################################################################################
# if len(sys.argv) < 2 or len(sys.argv) > 2:
	# print "Input argument error: The number of input arguments cannot be different from one!"
	# sys.exit(1)
# try:
	# f = open(str(sys.argv[1]),'r')
# except IOError:
	# print "Input data file error: The file '%s' was not found!" % str(sys.argv[1])
	# sys.exit(1)
# else:
	# f.close()
# sys.stdout.write("\r")
# sys.stdout.write("Processing: %4.1f%%" % 0.0)
# sys.stdout.flush()
topol_filename = str(sys.argv[1])
# #################################################################################################
# Initialize the parameters and the output variables
# #################################################################################################
# Initialize the basic output data

#number of days
if len(sys.argv) > 1:
        if '--n_days' in sys.argv:
            index_n_days=sys.argv.index('--n_days')
            n_days=int(sys.argv[index_n_days+1])
        else:
            n_days=1

n_simulations = int(3600 * 24)*n_days# Number of time samples (3600 sec. x 24 hours)

Hsmax = 0							# Initialize the number oh houses
arq = open(topol_filename,'r')
while True:
	C = arq.readline()
	if C == '':
		break
	while C.isspace() or C[0] == '*':
		C = arq.readline()
	if C == '':
		break
	C = [A.strip() for A in C.split(',')]
	if C[0] == 'HOUSE':
		Hsmax = Hsmax + 1
arq.close()
if Hsmax == 0:
	print "Input file error: No houses were found in the input data file!"
	sys.exit(1)
from electrical_data import general			# Establishes the harmonic input parameters
Hmax = (general['harm_max'] + general['harm_min'])/general['harm_step']
# Initial the house data array
from definitions import Zeros
VNeutral_Hs = [[]] * int(Hmax)
VPhase_Hs = [[]] * int(Hmax)
INeutral_Hs = [[]] * int(Hmax)
IPhase_Hs = [[]] * int(Hmax)
S_Hs = Zeros(2 * Hsmax, n_simulations)
# Q_Hs = Zeros(2 * Hsmax, n_simulations)
VNeutral_Hs = [Zeros(Hsmax, n_simulations) for A in VNeutral_Hs]
VPhase_Hs = [Zeros(2 * Hsmax, n_simulations) for A in VNeutral_Hs]
INeutral_Hs = [Zeros(Hsmax, n_simulations) for A in VNeutral_Hs]
IPhase_Hs = [Zeros(2 * Hsmax, n_simulations) for A in VNeutral_Hs]
# Initialize the transformer data array
VNeutral_TF = [[]] * int(Hmax)
VPhase_TF = [[]] * int(Hmax)
IPhase_TF = [[]] * int(Hmax)
S_TF = Zeros(2,n_simulations)
VNeutral_TF = [Zeros(1, n_simulations) for A in VNeutral_TF]
VPhase_TF = [Zeros(2, n_simulations) for A in VPhase_TF]
IPhase_TF = [Zeros(3, n_simulations) for A in IPhase_TF]
# #################################################################################################
# Read and format the input and default variables
# #################################################################################################
from ReadInput import ReadInput
app = ReadInput(topol_filename)
if app.err_state == 1:
	sys.exit(1)
# #################################################################################################
# Create the time-of-use and power profile of each appliance
# #################################################################################################
from definitions import CreateTimeOfUseProfile

app = CreateTimeOfUseProfile(app,n_simulations,n_days)
NodalVoltage = Zeros(len(app.NodesMap) + len(app.VS), int(Hmax))
for time_step in range(n_simulations):
	harm_order = app.harm_min
	old_state = [0.0] * (len(app.NodesMap) + len(app.VS))
	solution = 0
	it = [0] * int(Hmax)
	while harm_order <= app.harm_max:
		frequency = harm_order * 60.0
		it[solution] = 1
		# ########################################################################################
		# Power flow calculation
		# ########################################################################################
		from definitions import PrepareSolutionConditions
		from definitions import PowerFlow
		from definitions import PostPowerFlow
		while True:
			# Prepare the house equivalent impedance/current source for network power flow
			app = PrepareSolutionConditions(app, harm_order, NodalVoltage, old_state, time_step)
			# Solve the network power flow
			state, result_add = PowerFlow(app, frequency)
			# Update the house internal voltage for the next iteration
			if harm_order == 1:
				app = PostPowerFlow(app, state)
			# Stopping criteria
			Dif = [0.0] * len(state)
			for k in range(len(state)):
				Dif[k] = abs(state[k] - old_state[k])
			if max(Dif) < 0.05 or harm_order != 1 or it[solution] == 3:
				break
			old_state = state
			it[solution] = it[solution] + 1
		# ########################################################################################
		# Arrange the output data
		# ########################################################################################
		h = int((harm_order + 1.0) / 2.0) - 1
		for k in range(len(app.NodesMap) + len(app.VS)):
			NodalVoltage[k][solution] = state[k]
		for k in range(Hsmax):
			VNeutral_Hs[h][k][time_step] = result_add['VNeutral_Hs'][k]
			INeutral_Hs[h][k][time_step] = result_add['INeutral_Hs'][k]
		for k in range(2 * Hsmax):
			VPhase_Hs[h][k][time_step] = result_add['VPhase_Hs'][k]
			IPhase_Hs[h][k][time_step] = result_add['IPhase_Hs'][k]
			if harm_order == 1:
				S_Hs[k][time_step] = result_add['P_Hs'][k] + 1j*result_add['Q_Hs'][k]
				# Q_Hs[k][time_step] = result_add['Q_Hs'][k]
		VNeutral_TF[h][0][time_step] = state[4]
		VPhase_TF[h][0][time_step] = result_add['PhaseV_trafo'][0]
		VPhase_TF[h][1][time_step] = result_add['PhaseV_trafo'][1]
		IPhase_TF[h][0][time_step] = result_add['I_trafo'][0]
		IPhase_TF[h][1][time_step] = result_add['I_trafo'][1]
		IPhase_TF[h][2][time_step] = result_add['I_trafo'][2]
		if harm_order == 1:
			S_TF[0][time_step] = result_add['PQ_trafo'][0]
			S_TF[1][time_step] = result_add['PQ_trafo'][1]
		# ########################################################################################
		# Prepare the parameters for the next harmonic order
		# ########################################################################################
		solution = solution + 1
		harm_order = harm_order + app.harm_step
	# ############################################################################################
	# Update the process signal
	# ############################################################################################
	if (time_step + 1) % (n_simulations/1000) == 0:
		sys.stdout.write("\r")
		sys.stdout.write("Processing: %4.1f%%" % float(100.0 * float(time_step)/float(n_simulations)))
		sys.stdout.flush()
sys.stdout.write("\n")
sys.stdout.flush()
# #################################################################################################
# Create output csv files
# #################################################################################################
from definitions import OutputFileCSV
from os import path, makedirs
if not path.exists(path.dirname(path.abspath(__file__)) + r'\results'):
    makedirs(path.dirname(path.abspath(__file__)) + r'\results')
filename = 'PQ_Hs.csv'
OutputFileCSV(filename, S_Hs)
# filename = 'Q_Hs.csv'
# OutputFileCSV(filename, Q_Hs)
filename = 'PQ_Trafo.csv'
OutputFileCSV(filename, S_TF)
harm_order = 1
while harm_order <= app.harm_max:
	filename = 'IPhase_' + str(harm_order) + 'harm.csv'
	h = int((harm_order + 1.0) / 2.0) - 1
	OutputFileCSV(filename, IPhase_Hs[h])
	filename = 'VPhase_' + str(harm_order) + 'harm.csv'
	OutputFileCSV(filename, VPhase_Hs[h])
	filename = 'INeutral_' + str(harm_order) + 'harm.csv'
	OutputFileCSV(filename, INeutral_Hs[h])
	filename = 'VNeutral_' + str(harm_order) + 'harm.csv'
	OutputFileCSV(filename, VNeutral_Hs[h])
	filename = 'VNeutral_TF_' + str(harm_order) + 'harm.csv'
	OutputFileCSV(filename, VNeutral_TF[h])
	filename = 'VPhase_TF_' + str(harm_order) + 'harm.csv'
	OutputFileCSV(filename, VPhase_TF[h])
	filename = 'IPhase_TF_' + str(harm_order) + 'harm.csv'
	OutputFileCSV(filename, IPhase_TF[h])
	harm_order = harm_order + app.harm_step