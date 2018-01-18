import os
from numpy import real, imag

def PostPowerFlow(app, state):
	import math
	for House in app.houses:
		# total number of appliances connected on the current house
		appliance = len(House['appliances'])
		# room_number each appliance is connected to.
		Rmnum = [0] * appliance
		for Appliance in House['appliances']:
			# actualizes distances from one receptacle to the next one
			Rmnum[Appliance['app_num'] - 1] = Appliance['room_number'] # room_number has been identified by the previous subprograms
		# ####################### Iteration for single house impedance ####################
		for j in range(max(Rmnum)):
			# Loops over all rooms
			App_room = []
			# Determine the appliances connected in room j
			for Appliance in House['appliances']: 
				if Appliance['room_number'] - 1 == j:
					App_room.append(Appliance)
			conn = App_room[0]['app_conn'] - 1			# Room connection: A(=1), B(=2), AB(=3)
			ik = House['ifromBus'][conn] - 1
			im = House['itoBus'][conn] - 1
			# ################### Iteration for single room impedance ####################
			Vrcpl = state[ik] - state[im]			# The initial receptacle voltage is service panel voltage
			for Appliance in App_room: 
				Zaggr = Appliance['Zaggr']
				Zwr = Appliance['Zwr']
				Vrcpl = Vrcpl * Zaggr / (Zaggr + Zwr)
				Appliance['Vrcpl'] = Vrcpl			# Update the receptacle voltage
	return app

def PowerFlow(app, frequency):
	import math
	result_add = {}
	# Prepare the basic parameters for the load flow
	base_freq = 60.0
	tol = 1e-2
	Y = BuildYbus(app)
	index = len(Y)
	Ybus = Zeros((index + len(app.VS)), (index + len(app.VS)))
	#Ybus = [[0.0] * (index + len(app.VS))] * (index + len(app.VS))
	I = [0.0] * (index + len(app.VS))
	for i in range(index):
		for j in range(index):
			Ybus[i][j] = Y[i][j]
	# Include voltage sources in the circuit
	for VS in app.VS:
		Ybus[VS['ifromBus'][0] - 1][index] = 1.0
		Ybus[VS['itoBus'][0] - 1][index] = -1.0
		Ybus[index][VS['ifromBus'][0] - 1] = 1.0
		Ybus[index][VS['itoBus'][0] - 1] = -1.0
		if VS['freq'] > 0.0 and abs(VS['freq'] - frequency) > tol * base_freq:
			I[index] = 0.0
		else:
			I[index] = VS['mag'] * (math.cos(3.1416 / 180.0 * VS['ang']) + 1j * math.sin(3.1416 / 180.0 * VS['ang']))
			# I[index] = VS['mag'] * (math.cos(VS['ang']) + 1j * math.sin(VS['ang']))
		index = index + 1
	for House in app.houses:
		I[House['ifromBus'][0] - 1] = I[House['ifromBus'][0] - 1] - House['Ia']
		I[House['itoBus'][0] - 1] = I[House['itoBus'][0] - 1] + House['Ia']
		I[House['ifromBus'][1] - 1] = I[House['ifromBus'][1] - 1] - House['Ib']
		I[House['itoBus'][1] - 1] = I[House['itoBus'][1] - 1] + House['Ib']
		I[House['ifromBus'][2] - 1] = I[House['ifromBus'][2] - 1] - House['Iab']
		I[House['itoBus'][2] - 1] = I[House['itoBus'][2] - 1] + House['Iab']
	# Fill ground connection
	Ybus[0] = [0.0] * index
	Ybus[0][0] = 1.0
	I[0] = 0.0
	# Solve the nodal voltage equation
	# print Ybus
	# print I
	import numpy as np
	from numpy.linalg import inv, svd
	A = np.array(Ybus)
	B = np.array(I)
	# U,s,V = svd(A)
	# Ainv = np.dot(np.dot(V.T,inv(np.diag(s))),U.T)
	Ainv = inv(A)
	C = Ainv.dot(B)
	state = C.tolist()
	# print state
	horder = frequency/60.0
	# Calculate house voltage, current and power
	Hsnum = len(app.houses)
	PhaseV_Hs = [0.0] * (2 * Hsnum)
	PhaseI_Hs = [0.0] * (2 * Hsnum)
	NeutralV_Hs = [0.0] * Hsnum
	NeutralI_Hs = [0.0] * Hsnum
	Vphase = [0.0, 0.0]
	P_Hs = [0.0] * (2 * Hsnum)
	Q_Hs = [0.0] * (2 * Hsnum)
	for i in range(Hsnum):
		# Calculate the service panel current (phase A, phase B, and neutral)
		differ = 3 * i
		dV = [	state[7] - state[10 + differ], 
				state[8] - state[11 + differ],
				state[9] - state[12 + differ]]
		horder = frequency/60.0
		Z_Hs = app.branches[4 + 5 * i]['Z'] # Read the PCC-House impedance matrix
		Z_Hs_h = [[ Element.real + 1j * horder * Element.imag for Element in row] for row in Z_Hs]
		A = np.array(Z_Hs_h)
		B = np.array(dV)
		# U,s,V = svd(A)
		# Ainv = np.dot(np.dot(V.T,inv(np.diag(s))),U.T)
		Ainv = inv(A)
		C = Ainv.dot(B)
		Itemp = C.tolist() # Item = [Ia, In, Ib]
		PhaseI_Hs[2*i:2*i+2] = [Itemp[0], Itemp[2]]
		NeutralI_Hs[i] = Itemp[1]
		# Calculate the service panel phase voltage (phase AN, phase BN and neutral)
		ik = app.houses[i]['ifromBus'][0] - 1 # Room fromBus [11, 12, 11]...
		im = app.houses[i]['itoBus'][0] - 1	# Room toBus [12, 13, 13]...
		PhaseV_Hs[2*i] = state[ik] - state[im] # [V_an V_bn]
		ik = app.houses[i]['ifromBus'][1] - 1 # Room fromBus [11, 12, 11]...
		im = app.houses[i]['itoBus'][1] - 1	# Room toBus [12, 13, 13]...
		PhaseV_Hs[2*i+1] = state[ik] - state[im] # [V_an V_bn]
		NeutralV_Hs[i] = state[11+3*i]
		if horder == 1.0:
			Vphase = PhaseV_Hs[2*i:2*i+2]
			P_Hs[2*i] = Vphase[0] * np.conj(PhaseI_Hs[2*i])
			P_Hs[2*i] = P_Hs[2*i].real
			Q_Hs[2*i] = Vphase[0] * np.conj(PhaseI_Hs[2*i])
			Q_Hs[2*i] = Q_Hs[2*i].imag
			# The default Ib direction is from system to house, but the calculated Ib is the opposite, so we add a "-"
			P_Hs[2*i+1] = Vphase[1] * np.conj(PhaseI_Hs[2*i+1])
			P_Hs[2*i+1] = -P_Hs[2*i+1].real
			Q_Hs[2*i+1] = Vphase[1] * np.conj(PhaseI_Hs[2*i+1])
			Q_Hs[2*i+1] = -Q_Hs[2*i+1].imag
	result_add['VPhase_Hs'] = PhaseV_Hs
	result_add['IPhase_Hs'] = PhaseI_Hs#[Element.conjugate() for Element in PhaseI_Hs]
	result_add['VNeutral_Hs'] = NeutralV_Hs
	result_add['INeutral_Hs'] = NeutralI_Hs#[Element.conjugate() for Element in NeutralI_Hs]
	result_add['P_Hs'] = P_Hs
	result_add['Q_Hs'] = Q_Hs
	
	#print PhaseI_Hs
	#print NeutralI_Hs
	#print PhaseV_Hs
	#print NeutralV_Hs
	#print P_Hs
	#print Q_Hs
	
	# Calculate transformer voltage current, and power (secondary side)
	dV = [	state[5] - state[7],
			state[4] - state[8],
			state[6] - state[9]]
	Z_Hs = app.branches[3]['Z'] # Read the PCC-House impedance matrix
	Z_Hs_h = [[ Element.real + 1j * horder * Element.imag for Element in row] for row in Z_Hs]
	A = np.array(Z_Hs_h)
	B = np.array(dV)
	Ainv = inv(A)
	# U,s,V = svd(A)
	# Ainv = np.dot(np.dot(V.T,inv(np.diag(s))),U.T)
	C = Ainv.dot(B)
	Itemp = C.tolist() # Item = [Ia, In, Ib]
	I_trafo = Itemp
	# Calculate transformer power (secondary side)
	PQ_trafo = [0.0, 0.0]
	if horder == 1.0:
		PQ_trafo[0] = (state[5] - state[4]) * np.conj(Itemp[0])
		PQ_trafo[1] = (state[4] - state[6]) * np.conj(-Itemp[2])
	result_add['PhaseV_trafo'] = [state[5] - state[4], state[4] - state[6]]
	result_add['I_trafo'] = I_trafo
	result_add['PQ_trafo'] = PQ_trafo
	return state, result_add

def BuildYbus(app):
	totalbranches = 0
	tol = 1e-2
	# calculates number of branches that must be input in Ybus matrix (including nonzero shunts)
	for Branch in app.branches:
		totalbranches = totalbranches + Branch['phase']
		if Branch['phase'] == 1 and abs(Branch['Ysh'][0]) > tol:
			totalbranches = totalbranches + 2
		elif Branch['phase'] > 1 and max(max([[abs(Element) for Element in row] for row in Branch['Ysh']])) > tol:
			totalbranches = totalbranches + 2 * Branch['phase']
	A = Zeros(totalbranches, len(app.NodesMap))
	yprim = Zeros(totalbranches, totalbranches)
	#A = [[0.0] * len(app.NodesMap)] * totalbranches
	#yprim = [[0.0] * totalbranches] * totalbranches
	index = 0
	for Branch in app.branches:
		if Branch['phase'] == 1:
			A[index][Branch['ifromBus'][0] - 1] = 1.0
			A[index][Branch['itoBus'][0] - 1] = -1.0
			yprim[index][index] = Branch['Y'][0]
			index = index + 1
		elif Branch['phase'] > 1:
			for j in range(Branch['phase']):
				A[index][Branch['ifromBus'][j] - 1] = 1.0
				A[index][Branch['itoBus'][j] - 1] = -1.0
				for k in range(Branch['phase']):
					yprim[index][index + k - j] = Branch['Y'][j][k]
				index = index + 1
		# In case shunt elemnts are nonzero on PI branches
		if Branch['phase'] == 1 and abs(Branch['Ysh'][0]) > tol:
			A[index][Branch['ifromBus'][0] - 1] = 1.0
			A[index][0] = -1.0
			A[index + 1][Branch['itoBus'][0] - 1] = 1.0
			A[index + 1][0] = -1.0
			yprim[index][index] = Branch['Ysh'][0]
			yprim[index+1][index+1] = Branch['Ysh'][0]
			index = index + 2
		elif Branch['phase'] > 1 and max(max([[abs(Element) for Element in row] for row in Branch['Ysh']])) > tol:
			for j in range(Branch['phase']):
				A[index][Branch['ifromBus'][j] - 1] = 1.0
				A[index][0] = -1.0
				A[index + Branch['phase']][Branch['itoBus'][j] - 1] = 1.0
				A[index + Branch['phase']][0] = -1.0
				for k in range(Branch['phase']):
					yprim[index][index + k - j] = Branch['Ysh'][j][k]
					yprim[index + Branch['phase']][index + Branch['phase'] + k - j] = Branch['Ysh'][j][k]
				index = index + 1 + Branch['phase']
	#print yprim
	import numpy as np
	A = np.array(A)
	yprim = np.array(yprim)
	Ybus = yprim.dot(A)
	Ybus = A.T.dot(Ybus)
	Ybus = Ybus.tolist()
	#print Ybus
	#print app.NodesMap
	return Ybus

def PrepareSolutionConditions(app, harm_order, NodalVoltage, old_state, time_step):
	# This function actualizes house equivalent as well as other system parameters
	# for all frequancies, including fundamental
	import math
	Vinit = [120, 120, 240] # Initial guess for phase A, B, and AB voltages
	for House in app.houses:
		# phase and neutral cable parameters (ohm/length) inside the house
		Zphase = House['Zphase'].real + 1j * harm_order * House['Zphase'].imag
		Zneutral = House['Zneutral'].real + 1j * harm_order * House['Zphase'].imag
		# lower cable impedance: neutral parameters for AN and BN connections and phase parameters for AB connection 
		Zwr_pu = [Zphase+Zneutral, Zphase+Zneutral, 2*Zphase]
		# total number of appliances connected on the current house
		appliance = len(House['appliances'])
		# room_number each appliance is connected to.
		Rmnum = [0] * appliance
		dstn0 = [0] * appliance
		dstn = [0] * appliance
		for Appliance in House['appliances']:
			# actualizes distances from one receptacle to the next one
			Rmnum[Appliance['app_num'] - 1] = Appliance['room_number'] # room_number has been identified by the previous subprograms
			dstn0[Appliance['app_num'] - 1] = Appliance['dist']
			if Appliance['app_num'] == 1 or (Appliance['app_num'] > 1 and Rmnum[Appliance['app_num'] - 1] != Rmnum[Appliance['app_num'] - 2]):
				dstn[Appliance['app_num'] - 1] = dstn0[Appliance['app_num'] - 1]
			else:
				dstn[Appliance['app_num'] - 1] = dstn0[Appliance['app_num'] - 1] - dstn0[Appliance['app_num'] - 2]
		# house equivalent admittance for each phase (Za, Zb, Zab)
		Yeqq = [0.0] * 3
		# house equivalent harmonic current (Ia, Ib, Iab)
		Ieqq = [0.0] * 3
		# #############################################################
		# Iteration for single house impedance
		# #############################################################
		for j in range(max(Rmnum)):
			# Loops over all rooms
			Zaggr = 1e5 # room equivalent impedance
			Iaggr = 0.0 # room equivalent current
			App_room = []
			# Determine the appliances connected in room j
			for Appliance in House['appliances']: 
				if Appliance['room_number'] - 1 == j:
					App_room.append(Appliance)
			for k in range(len(App_room)):
				k_total = len(App_room) - 1 - k
				data_tmp = App_room[k_total]
				name = data_tmp['app_name']

				new_index=time_step%(3600*24/60)
				powratio = data_tmp['powratio'][new_index]

				switch_on = data_tmp['timeofuse'][time_step]
				if k_total < len(App_room) - 1:
					Zwr = Zwr_pu[App_room[0]['app_conn'] - 1] * dstn[App_room[k_total + 1]['app_num'] - 1]
				else:
					Zwr = 0
				if powratio < 0.01:
					powratio = 0.01
				# Part 1: Read the power data from app data package
				# If the appliance is on:
				if switch_on == 1.0:
					if name == 'RFR' or name == 'ASDFR':
						# the 3-sec spike is in progressing, the 2nd or more power flow iteration in this time step
						pfratio = data_tmp['pfratio'][time_step]
						Pld = powratio * data_tmp['app_qtt'] * data_tmp['P']
						Qld = Pld * math.sqrt(1 - pfratio**2)/pfratio
					elif name == 'PV':
						Pld = - powratio
						Qld = data_tmp['P'] * math.sqrt(1 - data_tmp['pf']**2)/data_tmp['pf']
					else:
						Pld = powratio * data_tmp['app_qtt'] * data_tmp['P']
						Qld = Pld * math.sqrt(1 - data_tmp['pf']**2)/data_tmp['pf']
					#print House['num_house'], data_tmp['room_number'], name
					# Part 2: Calculate the harmonic current injection
					if harm_order != 1 and data_tmp['loadcode'] != 0: # loadcode != 0 means nonliner load
						# print "epa!!!"
						Zaggr = Zaggr + Zwr
						index = int((harm_order - app.harm_min)/2.0)# Harmonic under under study
						Vrcpl = data_tmp['Vrcpl'] # Receptacle voltage
						Iload = (Pld + 1j * Qld) / Vrcpl
						Iload = Iload.conjugate()
						if data_tmp['loadcode'] == 1:
							# Using default spectrum
							Iapp = app.spectrum[data_tmp['app_name']]
						elif data_tmp['loadcode'] == -1:
							# Using external spectrum
							Iapp = data_tmp['external_espectrum']
						# Spectrum data for the current appliance
						Ispectr_mag = Iapp[0][index]
						Ispectr_ang = Iapp[1][index]
						if name == 'PV':
							Irated = - data_tmp['P'] / 240.0
							Ih_mag = abs(Irated) * Ispectr_mag / Iapp[0][0]
							Ih_ang = Ispectr_ang * 3.1416 / 180.0 + harm_order * (math.atan2(Iload.imag, Iload.real) - Iapp[1][0] * 3.1416 / 180.0)
							Ihload = - Ih_mag * (math.cos(Ih_ang) + 1j * math.sin(Ih_ang))
							Iaggr = Iaggr + Ihload
							# Calculate the Northon equivalent impedance
							Z1 = 0.05 + 1j * harm_order * 377 * 2e-3
							Zc = 2.0 + 1.0 / (1j * harm_order * 377 * 10.4e-6)
							Z2 = 0.05 + 1j * harm_order * 377 * 1.3e-3
							Zinv = Z1 * Zc / (Z1 + Zc)
							# Integrate the Norton impedance to the PCC
							if k_total == len(App_room) - 1:
								# If there is already load integrated before, this load is a parallel impedance
								Zaggr = Zinv + Z2
							else:
								Zaggr = (Zaggr+Zwr)*Zinv/((Zaggr+Zwr)+Zinv)+Z2
						else:
							Ih_mag = math.sqrt(Iload.real**2 + Iload.imag**2) * Ispectr_mag / Iapp[0][0]
							Ih_ang = Ispectr_ang * 3.1416 / 180.0 + harm_order * (math.atan2(Iload.imag, Iload.real) - Iapp[1][0] * 3.1416 / 180.0)
							Ihload = - Ih_mag * (math.cos(Ih_ang) + 1j * math.sin(Ih_ang))
							Iaggr = Iaggr + Ihload
					else: 
						# In case harm_order==1 (fundamental freq.) or loadcode==0 (linear code)
						# calculates equivalent impedance for fundamental frequency and linear loads on harmonic frequencies.
						sum1 = 0.0
						sum2 = 0.0
						sum3 = 0.0
						for Element in NodalVoltage:
							sum1 = sum1 + math.sqrt(Element[0].real**2 + Element[0].imag**2)
						for Element in old_state:
							sum2 = sum2 + math.sqrt(Element.real**2 + Element.imag**2)
						for Element1 in NodalVoltage:
							for Element2 in Element1:
									sum3 = sum3 + math.sqrt(Element2.real**2 + Element2.imag**2)
						if sum1 < 0.1 and sum2 < 0.1:
							# First iteration of fundamental load flow
							Rld = Vinit[App_room[0]['app_conn'] - 1]**2 * Pld / (Pld**2 + Qld**2)
							Xld = harm_order * Vinit[App_room[0]['app_conn'] - 1]**2 * Qld / (Pld**2 + Qld**2)
							# print "Tensoes 1:"
							# print Vinit
						elif sum3 < 0.1:
							# subsequent iterations of fundamental load flow
							Vrcpl = data_tmp['Vrcpl']
							Rld = (Vrcpl.real**2 + Vrcpl.imag**2) * Pld / (Pld**2 + Qld**2)
							Xld = harm_order * (Vrcpl.real**2 + Vrcpl.imag**2) * Qld / (Pld**2 + Qld**2)
							# print "Tensoes 1:"
							# print Vinit
						else:
							# harmonic load flows
							Vrcpl = data_tmp['Vrcpl']
							Rld = (Vrcpl.real**2 + Vrcpl.imag**2) * Pld / (Pld**2 + Qld**2)
							Xld = harm_order * (Vrcpl.real**2 + Vrcpl.imag**2) * Qld / (Pld**2 + Qld**2)
						Zld = Rld + 1j * Xld
						if k_total == len(App_room) - 1:
							Zaggr = Zld
						else:
							Zaggr = (Zaggr + Zwr) * Zld / ((Zaggr + Zwr) + Zld)
				# If the appliance is off:
				else:
					# if the appliance is not ON: it is open-circuit and will not be processed.
					Zaggr = Zaggr + Zwr
					flag_on = 0
				App_room[k]['Zwr'] = Zwr
				App_room[k]['Zaggr'] = Zaggr
			# actualize impedance and current values considering different rooms in parallel
			if Zaggr != 0.0:
				Yeqq[App_room[0]['app_conn'] - 1] = 1 / (Zaggr + Zwr_pu[App_room[0]['app_conn'] - 1] * dstn[App_room[0]['app_num'] - 1]) + Yeqq[App_room[0]['app_conn'] - 1]
			Ieqq[App_room[0]['app_conn'] - 1] = Ieqq[App_room[0]['app_conn'] - 1] + Iaggr
		# #############################################################
		# store equivalent parameters for the correspondent house
		# #############################################################
		House['Ia'] = Ieqq[0]
		House['Ib'] = Ieqq[1]
		House['Iab'] = Ieqq[2]
		app.branches[House['ID'] - 1]['Y'][0] = Yeqq[0]
		app.branches[House['ID'] + 0]['Y'][0] = Yeqq[1]
		app.branches[House['ID'] + 1]['Y'][0] = Yeqq[2]
	# primary system equivalent impedance must be treated differently
	Zld = [[app.equivalent['Zpp'][harm_order - 1], app.equivalent['Zpn'][harm_order - 1]], [app.equivalent['Zpn'][harm_order - 1], app.equivalent['Znn'][harm_order - 1]]]
	import numpy as np
	from numpy.linalg import inv
	A = np.array(Zld)
	Ainv = inv(A)
	app.branches[0]['Y'] = Ainv.tolist()
	app.branches[0]['Ysh'] = [[element / harm_order for element in row] for row in app.branches[0]['Ysh']]
	for i in range(1,len(app.branches)):
		if app.branches[i]['type'] == 'HOUSE_Z':
			continue
		elif app.branches[i]['type'] == 'PI':
			if app.branches[i]['phase'] == 3:
				Zld = [[element.real + 1j * harm_order * element.imag for element in row] for row in app.branches[i]['Z']]
				A = np.array(Zld)
				Ainv = inv(A)
				app.branches[i]['Y'] = Ainv.tolist()
				app.branches[i]['Ysh'] = [[element / harm_order for element in row] for row in app.branches[i]['Ysh']]
			elif app.branches[i]['phase'] == 1:
				Zld = app.branches[i]['Z'][0].real + 1j * harm_order * app.branches[i]['Z'][0].imag
				app.branches[i]['Y'][0] = 1.0 / Zld
				app.branches[i]['Ysh'][0] = app.branches[i]['Ysh'][0] / harm_order
		elif app.branches[i]['type'] == 'TRAFO':
			app.branches[i]['Y'] = [[element / harm_order for element in row] for row in app.branches[i]['Yprim']]
			app.branches[i]['Ysh'] = [[element / harm_order for element in row] for row in app.branches[i]['Ysh']]
	#for Branch in app.branches:
	#	print Branch['Y']
	#for House in app.houses:
	#	print House['Ia'], House['Ib'], House['Iab']
	return app

def instant_generation(begin_time,end_time):
    """
        return a instant between begin_time and end_time
        ----
        example: 
        
        input:
        
        begin_time='6:00'
        end_time='7:30'
        
        output:
        
        421  #--- 421 is equal '7:01'
        
    """ 
    [hour_begin,minutes_begin]=str.split(begin_time,':')
    [hour_end,minutes_end]=str.split(end_time,':')
    
import datetime
import random
import math
import matplotlib.pyplot as plt

def instant_generation(begin_time,end_time):
    """
        return a instant between begin_time and end_time
        ----
        example: 
        
        input:
        
        begin_time='6:00'
        end_time='7:30'
        
        output:
        
        421  #--- 421 is equal '7:01'
        
    """
    midnight=datetime.datetime.strptime('00:00', "%H:%M")
    
    b = datetime.datetime.strptime(begin_time, "%H:%M") 
    e = datetime.datetime.strptime(end_time, "%H:%M")
    
    return random.randint(delta_in_minutes(b,midnight),delta_in_minutes(e,midnight))

def delta_in_minutes(t1,t2):
    return math.ceil((t1-t2).total_seconds()/60)


def generate_occupancy(Occupancy_pro,interval_begin,interval_finish,n_days=1,correct=False):
    
    t_begin = instant_generation(interval_begin[0],interval_begin[1])
    t_finish = instant_generation(interval_finish[0],interval_finish[1])
    
    if correct:
         if t_finish < t_begin + 20:
            t_finish = t_begin + 20
    
    seconds_in_day=24*60
    for i in range(n_days):
        Occupancy_pro[t_begin+i*seconds_in_day - 1:t_finish+i*seconds_in_day - 1] = [1.0] * (t_finish - t_begin)
    
    return Occupancy_pro
    

def deal_with_work_type(Occupancy_pro,work_type,type_house_dictionary,n_days):
    
    dictionary=type_house_dictionary[work_type]
    
    keys=dictionary.keys()
    
    if ('waking' not in keys) & ('leaving' in keys) &  ('arriving' in keys) & ('sleeping' not in keys):
        
        Occupancy_pro=generate_occupancy(Occupancy_pro,dictionary['arriving'],dictionary['leaving'],n_days)
    
    elif ('waking' in keys) & ('leaving' not in keys) &  ('arriving' not in keys) & ('sleeping' in keys):
        
        Occupancy_pro=generate_occupancy(Occupancy_pro,dictionary['waking'],dictionary['sleeping'],n_days)
    
    else:
        
        Occupancy_pro=generate_occupancy(Occupancy_pro,dictionary['waking'],dictionary['leaving'],n_days,correct=True)
        Occupancy_pro=generate_occupancy(Occupancy_pro,dictionary['arriving'],dictionary['sleeping'],n_days)
      
    return Occupancy_pro
    
def start_end(startpoint,m,cycle_duration,min_duration):
    
    starting = int(60 * (startpoint + (m * cycle_duration)) + math.floor(random.random() * 60))
    ending  = int(60 * (startpoint + min_duration + (m * cycle_duration)) + math.floor(random.random() * 60))
    
    return starting,ending

def add_dynamic(duration_dynamic,Appliance,name,state):
    
    Dyn = Dynamic_check(Appliance, name)
    for n_dyn in range(duration_dynamic):
        powratio[starting + n_dyn] = Dyn['Pstrt'][n_dyn]
        states[starting + n_dyn] = [state]*len(Dyn['Pstrt'][n_dyn])
        pfratio[starting + n_dyn] = Dyn['PFstrt'][n_dyn]
    
    
from TimeOfUseProfile import house_work_types as house_wt

def CreateTimeOfUseProfile(app,n_simulations,n_days):
    
	import math
	import random    
	defined_number=n_simulations/60
    
	day_type = app.day_type
	Irr_full = app.Irr_full
	Irr_pu = PVsimulation(day_type, Irr_full,defined_number)
	
	    
	for House in app.houses:
		Occupancy_pro = [0.0] * defined_number
		k = House['size']/2.5 # Ojo!! Numero magico (Number of people)
		
		Occupancy_pro=deal_with_work_type(Occupancy_pro,House['work_type'],house_wt.work_type,n_days)
		
		# WSH and DRY, PC and LCD are correlated, so getting them appart from the other 20+ independent appliances
		ind_wsh = []
		ind_dry = []
		ind_pc = []
		ind_lcd = []
		for Appliance in House['appliances']:
			Appliance['processed'] = False
			if Appliance['app_name'] == 'PC':
				ind_pc.append(Appliance['app_num'])
			elif Appliance['app_name'] == 'LCD':
				ind_lcd.append(Appliance['app_num'])
			elif Appliance['app_name'] == 'WSH':
				ind_wsh.append(Appliance['app_num'])
			elif Appliance['app_name'] == 'DRY' or Appliance['app_name'] == 'ASDDRY':
				ind_dry.append(Appliance['app_num'])
		enter_wshdry = 0
		enter_pc = 0
		P_noise = House['noise'];

		for Appliance in House['appliances']:
			name = Appliance['app_name']
			if not name == 'PV' and not name == 'PHEV':
				# PV and PHEV have nothing to do with operation cycles
				n_cycles = int(math.ceil(app.appliance_char[name]['num_switchONs']))
				min_duration = app.appliance_char[name]['min_cycleduration']
				max_duration = app.appliance_char[name]['max_cycleduration']
			# Initialize time-of-use profile
			timesec0 = [0] * n_simulations#(3600 * 27)
			timesec = [0] * n_simulations
			powratio0 = [0.01] * n_simulations#(3600 * 27)
			powratio = [0.01] * n_simulations
			states0 = [0.01] *  n_simulations#(3600 * 27)
			states = [0.01] * n_simulations
			timesec_wsh0 = [0]  * n_simulations#(3600 * 27)
			timesec_wsh = [0] * n_simulations
			powratio_wsh0 = [0.01] * n_simulations#(3600 * 27)
			powratio_wsh = [0.01] * n_simulations
			states_wsh0 = [0.01] * n_simulations#(3600 * 27)
			states_wsh = [0.01] * n_simulations            
			pfratio = [0] * n_simulations
			Power_PV_sec =[0] * n_simulations
			time_PHEV_sec = [0] * n_simulations
			# Switch-on events are assigned probabilistically in seconds
            
            
            
			if Appliance['timeofusecode'] == 1:
				# ############################# Freezer ##########################
				if name == 'FRZR':
					cycle_duration = defined_number/n_cycles
					startpoint = int(math.floor(random.random() * (cycle_duration - min_duration)))
					for m in range(n_cycles):
						# Default: has a 6min cycle every 20 min interval during all day (6 min ON - 14 min OFF)
						starting 	= int(60 * (startpoint + (m * cycle_duration)) + math.floor(random.random() * 60))
						ending 		= int(60 * (startpoint + min_duration + (m * cycle_duration)) + math.floor(random.random() * 60))
						timesec[starting:ending] = [1] * (ending - starting)
						# Add spike
						Dyn = Dynamic_check(Appliance, name)
						for n_dyn in range(3):
							# When the fidge start each time, insert the 3 sec. dynamic process
							powratio[starting + n_dyn] = Dyn['Pstrt'][n_dyn]
							states[starting + n_dyn] = [1]*len(Dyn['Pstrt'][n_dyn])
							pfratio[starting + n_dyn] = Dyn['PFstrt'][n_dyn]
						# Add multi-stage
						T1 = round(random.normalvariate(30, 5))
						P1_ratio = round(random.normalvariate(1.2, 0.05))
						T1 = int(max([T1, 10]))
						powratio[starting + 3:starting + T1] = [P1_ratio] * (T1 - 3)
						powratio[starting + T1:ending] = [1.0] * (ending - starting - T1)
						states[starting + 3:starting + T1] = [1] * (T1 - 3)
						states[starting + T1:ending] = [2] * (ending - starting - T1)                        
						# Add fluctuation
						n_fluc = ending - starting - 90
						for n in range(n_fluc):
							num_fluc = math.floor(n/30)
							if (n % 60) < 30:
								powratio[starting + 60 + n] = 0.97 + (n - 30 * num_fluc) * 0.06 / 30
								states[starting + 60 + n] = 3                                
							else:
								powratio[starting + 60 + n] = 1.03 - (n - 30 * num_fluc) * 0.06 / 30
								states[starting + 60 + n] = 4                                 
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states                    
					Appliance['processed'] = True
					continue
				# ############################# Regular fridge or ASD fridge ##########################
				if name == 'RFR' or name == 'ASDFR':
					cycle_duration = defined_number/n_cycles
					startpoint = int(math.floor(random.random() * (cycle_duration - min_duration)))
					for m in range(n_cycles):
						# Default: has a 9 min cycle every 30 min interval during all day (9 min ON - 21 min OFF)
						starting 	= int(60 * (startpoint + (m * cycle_duration)) + math.floor(random.random() * 60))
						ending 		= int(60 * (startpoint + min_duration + (m * cycle_duration)) + math.floor(random.random() * 60))
						timesec[starting:ending] = [1] * (ending - starting)
						powratio[starting:ending] = [1] * (ending - starting)
						states[starting:ending] = [1] * (ending - starting)                        
						pfratio[starting:ending] = [Appliance['pf']] * (ending - starting)
						Dyn = Dynamic_check(Appliance, name)
						for n_dyn in range(3):
							# When the fidge start each time, insert the 3 sec. dynamic process
							powratio[starting + n_dyn] = Dyn['Pstrt'][n_dyn]
							pfratio[starting + n_dyn] = Dyn['PFstrt'][n_dyn]
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states                    
					Appliance['pfratio'] = pfratio
					Appliance['processed'] = True
					continue
				# ############################# Furnace ####################################
				if name == 'FUR':
					cycle_duration = defined_number/n_cycles
					startpoint = int(math.floor(random.random() * (cycle_duration - min_duration)))
					for m in range(n_cycles):
						# Number of subcycles can be 3, 4, or 5
						n_subcycles = 3 + int(math.floor(3 * random.random()))
						starting = int(60 * (startpoint + (m * cycle_duration)) + math.floor(random.random() * 60))
						timesec0[max([starting - 60, 0]):starting] = [1] * (starting - max([starting - 60, 0]))
						# the pre-starting process of fornance
						powratio0[max([starting - 60, 0]):starting] = [random.normalvariate(0.2, 0.01)] * (starting - max([starting - 60, 0]))
						states0[max([starting - 60, 0]):starting] = [1] * (starting - max([starting - 60, 0]))                        
						# powratio0[starting-2:starting] = [0.02] * 2
						for n in range(n_subcycles):
							T = int(round(random.normalvariate(360, 30)))
							Ton = int(round(T * (0.5 + 0.2 * random.random()))) # Duration of the motor-on
							timesec0[starting:starting + T] = [1] * T
							Phea_ratio = random.normalvariate(0.8, 0.01) # Power ratio of motor-off
							powratio0[starting:starting + Ton] = [1] * Ton
							powratio0[starting + Ton:starting + T - 2] = [Phea_ratio] * (T - 2 - Ton)
							powratio0[starting + T - 2:starting + T] = [1 - Phea_ratio] * 2
							states0[starting:starting + Ton] = [2] * Ton
							states0[starting + Ton:starting + T - 2] = [3] * (T - 2 - Ton)
							states0[starting + T - 2:starting + T] = [4] * 2                            
							pfratio[starting:starting + Ton] = [Appliance['pf']] * Ton
							pfratio[starting + Ton:starting + T] = [0.98] * (T - Ton) # Ojo, numero magico!!!
							# Add spike
							Dyn = Dynamic_check(Appliance, name)
							for n_dyn in range(3):
								# When the fidge start each time, insert the 3 sec. dynamic process
								powratio0[starting + n_dyn] = Dyn['Pstrt'][n_dyn]
								pfratio[starting + n_dyn] = Dyn['PFstrt'][n_dyn]
							starting = starting + T
					ending = starting
					timesec, powratio, states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states                    
					Appliance['pfratio'] = pfratio
					Appliance['processed'] = True
					continue
				prob_curve = Appliance['prob_curve']
				C = [0] * defined_number
				for i in range(defined_number):
					new_index=i%(3600*24/60)                    
					try:
						C[i] = prob_curve[new_index] * Occupancy_pro[i]
					except:
						print(i)
						print(len(prob_curve))
						print(len(Occupancy_pro))
						print('-*-')                        
				c = 100 / sum(C)
				# ############################# PC or LCD ####################################
				if name == 'PC' or name == 'LCD':
					if Appliance['processed'] == True:
						continue
					if len(ind_pc) != len(ind_lcd) and enter_pc + 1 > min([len(ind_pc), len(ind_lcd)]):
						ending = 0
						for time_step in range(defined_number):
							if Occupancy_pro[time_step] == 1.0 and timesec[60 * time_step] == 0:
								new_index=time_step%(24*60)                                 
								threshold = prob_curve[new_index] * n_cycles * k * c
								if random.random() * 100 < threshold:
									duration = min_duration + math.floor((max_duration - min_duration) * random.random())
									starting = int(60 * time_step + math.floor(random.random() * 60))
									ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
									timesec0[starting:ending] = [1.0] * (ending - starting)
									powratio0[starting:ending] = [1.0] * (ending - starting)
									states0[starting:ending] = [1] * (ending - starting)                                    
									# Add multi-stage
									T1 = int(round(60 * duration * (0.2 + 0.4 * random.random())))
									T2 = int(round(60 * duration * (0.2 + 0.5 * random.random())))
									P2_ratio = random.normalvariate(1.3, 0.03)
									powratio0[starting + T1:ending + T1 + T2] = [P2_ratio] * (ending + T2 - starting)
									powratio0[starting + T1:ending + T1 + T2] = [P2_ratio] * (ending + T2 - starting)
									state0[starting + T1:ending + T1 + T2] = [2] * (ending + T2 - starting)
									state0[starting + T1:ending + T1 + T2] = [2] * (ending + T2 - starting)                         
									for i in range(starting,ending):
										powratio0[i] = powratio0[i] + 0.04 * P_noise[i - starting]
						timesec, powratio, states = TimeModify2(timesec0, powratio0, states, ending)
						Appliance['timeofuse'] = timesec
						Appliance['powratio'] = powratio
						Appliance['states'] = states                        
						Appliance['pfratio'] = pfratio
						Appliance['processed'] = True
						enter_pc = enter_pc + 1
					else:
						ending = 0
						for time_step in range(defined_number):
							if Occupancy_pro[time_step] == 1.0 and timesec[60 * time_step] == 0:
								new_index=time_step%(24*60)                                  
								threshold = prob_curve[new_index] * n_cycles * k * c
								if random.random() * 100 < threshold:
									# PC
									duration = min_duration + math.floor((max_duration - min_duration) * random.random())
									starting = int(60 * time_step + math.floor(random.random() * 60))
									ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
									timesec0[starting:ending] = [1.0] * (ending - starting)
									powratio0[starting:ending] = [1.0] * (ending - starting)
									states0[starting:ending] = [1] * (ending - starting)                                    
									# Add multi-stage
									T1 = int(round(duration * (0.2 + 0.4 * random.random())))
									T2 = int(round(duration * (0.2 + 0.5 * random.random())))
									P2_ratio = random.normalvariate(1.3, 0.03)
									powratio0[starting + T1:ending + T1 + T2] = [P2_ratio] * (ending + T2 - starting)
									#powratio0[starting + T1:ending + T1 + T2] = [P2_ratio] * (ending + T2 - starting)
									states0[starting + T1:ending + T1 + T2] = [2] * (ending + T2 - starting)
									for i in range(starting,ending):
										powratio0[i] = powratio0[i] + 0.08 * P_noise[i - starting]
						timesec, powratio, states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
						House['appliances'][ind_pc[enter_pc] - 1]['timeofuse'] = timesec
						House['appliances'][ind_pc[enter_pc] - 1]['powratio'] = powratio
						House['appliances'][ind_pc[enter_pc] - 1]['states'] = states                        
						House['appliances'][ind_pc[enter_pc] - 1]['processed'] = True
						House['appliances'][ind_lcd[enter_pc] - 1]['timeofuse'] = timesec
						House['appliances'][ind_lcd[enter_pc] - 1]['powratio'] = powratio
						House['appliances'][ind_lcd[enter_pc] - 1]['states'] = states                        
						House['appliances'][ind_lcd[enter_pc] - 1]['processed'] = True
						enter_pc = enter_pc + 1
					continue
				# ################ Washer, dryer and ASD dryer ##########################
				if name == 'WSH' or name == 'DRY' or name == 'ASDDRY':
					if Appliance['processed'] == True:
						continue
					min_duration_wsh = app.appliance_char[House['appliances'][ind_wsh[enter_wshdry] - 1]['app_name']]['min_cycleduration']
					max_duration_wsh = app.appliance_char[House['appliances'][ind_wsh[enter_wshdry] - 1]['app_name']]['max_cycleduration']
					min_duration_dry = app.appliance_char[House['appliances'][ind_dry[enter_wshdry] - 1]['app_name']]['min_cycleduration']
					max_duration_dry = app.appliance_char[House['appliances'][ind_dry[enter_wshdry] - 1]['app_name']]['max_cycleduration']
					ending_wsh = 0
					ending_dry = 0
					for time_step in range(defined_number):
						if Occupancy_pro[time_step] == 1.0 and timesec[60 * time_step] == 0:
							new_index=time_step%(24*60)                            
							threshold = prob_curve[new_index] * n_cycles * k * c
							if random.random() * 100 < threshold:
								# Washing machine
								duration = min_duration + math.ceil((max_duration_wsh - min_duration_wsh) * random.random())
								starting_wsh = int(60 * time_step + math.floor(random.random() * 60))
								starting_wsh0 = starting_wsh
								ending_wsh = int(60 * (time_step + duration) + math.floor(random.random() * 60))
								timesec_wsh0[starting_wsh:ending_wsh] = [1.0] * (ending_wsh - starting_wsh)
								# Add pulse
								while True:
									T = int(round(random.normalvariate(30,1)))
									Ton = int(T * round(0.5 + 0.1 * random.random())) # The duration of motor-on
									powratio_wsh0[starting_wsh:starting_wsh + Ton] = [1.0] * Ton
									powratio_wsh0[starting_wsh + Ton:starting_wsh + Ton + T] = [0.01] * Ton
									states_wsh0[starting_wsh:starting_wsh + Ton] = [1] * Ton
									states_wsh0[starting_wsh + Ton:starting_wsh + Ton + T] = [2] * Ton                               
									starting_wsh = starting_wsh + T
									if starting_wsh > ending_wsh:
											break
								# Add noise
								for i in range(starting_wsh0, ending_wsh):
									powratio_wsh0[i] = powratio_wsh0[i] + 0.06 * P_noise[i - starting_wsh0]
								# dryer
								duration = int(math.floor(min_duration_dry + (max_duration_dry - min_duration_dry) * random.random()))
								starting_dry = int(60 * (time_step + 40) + math.floor(60 * random.random()))
								starting_dry0 = starting_dry
								ending_dry = int(60 * (time_step + duration + 40) + math.floor(random.random() * 60))
								timesec0[starting_dry:ending_dry] = [1.0] * (ending_dry - starting_dry)
								for n in range(starting_dry, ending_dry):
									powratio0[n] = powratio0[n] + (-0.05 + 0.06 * random.random())
									states0[n] = states0[n] + (1)                                    
								while True:
									T = int(round(random.normalvariate(240,20)))
									Ton = int(T * round(0.5 + 0.1 * random.random())) # The duration of motor-on
									for n_decr in range(10):
										powratio0[starting_dry + n_decr] = 1.07 - 0.007 * (n_decr - 1)
										states0[starting_dry + n_decr] = 2
									powratio0[starting_dry + 10:starting_dry + Ton] = [1.0] * (Ton - 10)
									powratio0[starting_dry + Ton:starting_dry + T] = [0.001] * (T - Ton)
									states0[starting_dry + 10:starting_dry + Ton] = [3] * (Ton - 10)
									states0[starting_dry + Ton:starting_dry + T] = [3] * (T - Ton)                                   
									starting_dry = starting_dry + T
									if starting_dry > ending_dry:
											break
								for i in range(starting_dry0, ending_dry):
									powratio0[i] = powratio0[i] + 0.005 * P_noise[i - starting_dry0]
					# Check if the operation time goes to the next day
					timesec_wsh, powratio_wsh,states_wsh = TimeModify2(timesec_wsh0, powratio_wsh0,states_wsh0, ending_wsh)
					# Check if the operation time goes to the next day
					timesec, powratio,states = TimeModify2(timesec0, powratio0, states0,ending_dry)
					# Add duty cycle to the DRY, 1 min per small cycle
					House['appliances'][ind_wsh[enter_wshdry] - 1]['timeofuse'] = timesec_wsh
					House['appliances'][ind_wsh[enter_wshdry] - 1]['powratio'] = powratio_wsh
					House['appliances'][ind_wsh[enter_wshdry] - 1]['states'] = states_wsh
					House['appliances'][ind_wsh[enter_wshdry] - 1]['processed'] = True
					House['appliances'][ind_dry[enter_wshdry] - 1]['timeofuse'] = timesec
					House['appliances'][ind_dry[enter_wshdry] - 1]['powratio'] = powratio
					House['appliances'][ind_dry[enter_wshdry] - 1]['states'] = states
					House['appliances'][ind_dry[enter_wshdry] - 1]['processed'] = True
					enter_wshdry = enter_wshdry + 1
					continue
				# ########################## Stove ################################
				if name == 'STO':
					ending = 0
					for time_step in range(defined_number):
						if Occupancy_pro[time_step] == 1.0 and timesec[60 * time_step] == 0:
							new_index=time_step%(3600*24/60)
							threshold = prob_curve[new_index] * n_cycles*n_days * k * c
							if random.random() * 100 < threshold:
								duration = min_duration + math.ceil((max_duration - min_duration) * random.random())
								starting = int(60 * time_step + math.floor(random.random() * 60))
								ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
								timesec0[starting:ending] = [1.0] * (ending - starting)
								Ton_ratio = 0.15 + 0.25 * random.random()		# The duration of motor-on
								while True:
									T = int(round(random.normalvariate(30,2)))
									Ton = int(round(Ton_ratio * T))
									powratio0[starting:starting + Ton] = [1.0] * Ton
									states0[starting:starting + Ton] = [1.0] * Ton
									powratio0[starting + Ton:starting + T] = [0.001] * (T - Ton)
									states0[starting + Ton:starting + T] = [2.0] * (T - Ton)
									starting = starting + T
									if starting > ending:
										break
					# Check if the operation time goes to the next
					timesec, powratio,states = TimeModify2(timesec0, powratio0,states0, ending,n_simulations)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states']=states
					Appliance['processed'] = True
					continue
				# ########################## Coffee maker, heater ##########################
				if name == 'COF' or name == 'HEA':
					ending = 0
					for time_step in range(defined_number):
						if Occupancy_pro[time_step] == 1.0 and timesec[60 * time_step] == 0:
							new_index=time_step%(3600*24/60)
							threshold = prob_curve[new_index] * n_cycles * k * c
							if random.random() * 100 < threshold:
								duration = min_duration + math.ceil((max_duration - min_duration) * random.random())
								starting = int(60 * time_step + math.floor(random.random() * 60))
								ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
								timesec0[starting:ending] = [1.0] * (ending - starting)
								if name == 'COF':
									T = int(round(random.normalvariate(30,2)))
									Ton = int(round(T * (0.3 + 0.4 * random.random())))
									while True:
										powratio0[starting:starting + Ton] = [1.0] * Ton
										powratio0[starting + Ton:starting + T] = [0.001] * (T - Ton)
										states0[starting:starting + Ton] = [1] * Ton
										states0[starting + Ton:starting + T] = [2] * (T - Ton)
										starting = starting + T
										if starting > ending:
											break
								T = int(round(random.normalvariate(40,2)))
								Ton = int(round(T * (0.1 + 0.4 * random.random())))
								while True:
									powratio0[starting:starting + Ton] = [1.0] * Ton
									powratio0[starting + Ton:starting + T] = [0.01] * (T - Ton)
									states0[starting:starting + Ton] = [1] * Ton
									states0[starting + Ton:starting + T] = [3] * (T - Ton)
									starting = starting + T
									if starting > ending:
										break
					timesec, powratio, states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states
					Appliance['processed'] = True
					continue
				# ########################## Vacuum Cleaner ##########################
				if name == 'VAC':
					for time_step in range(defined_number):
						if Occupancy_pro[time_step] == 1.0 and timesec[60 * time_step] == 0:
							new_index=time_step%(24*60)
							threshold = prob_curve[new_index] * n_cycles * k * c
							if random.random() * 100 < threshold:
								duration = min_duration + math.ceil((max_duration - min_duration) * random.random())
								starting = int(60 * time_step + math.floor(random.random() * 60))
								ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
								timesec0[starting:ending] = [1.0] * (ending - starting)
								powratio0 = timesec0
								powratio0[starting] = 1.5 - 0.25 * 2
								states0 = timesec0
								states0[starting] = 2
					timesec, powratio, states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states
					Appliance['processed'] = True
					continue
				# ########################## Laptop ##########################
				if name == 'LAP':
					ending = 0
					for time_step in range(defined_number):
						if Occupancy_pro[time_step] == 1.0 and timesec[60 * time_step] == 0 and timesec[60 * (time_step - 2)] == 0:
							new_index=time_step%(24*60)
							threshold = prob_curve[new_index] * n_cycles * k * c
							if random.random() * 100 < threshold:
								duration = min_duration + math.floor((max_duration - min_duration) * random.random())
								starting = int(60 * time_step + math.floor(random.random() * 60))
								ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
								timesec0[starting:ending] = [1.0] * (ending - starting)
								# Add multi-stage
								T1 = max([int(round(random.normalvariate(30,5))), 10])
								P1_ratio = max([random.normalvariate(1.4,0.05), 1.2])
								powratio0[starting:starting + T1] = [P1_ratio] * T1
								powratio0[starting + T1:ending] = [1.0] * (ending - starting - T1)
								states0[starting:starting + T1] = [1] * T1
								states0[starting + T1:ending] = [2] * (ending - starting - T1)
								T2 = int(round(60 * duration * (0.1 + 0.4 * random.random())))
								T3 = int(round(60 * duration * (0.2 + 0.5 * random.random())))
								P3_ratio = random.normalvariate(1.3,0.05)
								powratio0[starting + T1 + T2:starting + T1 + T2 + T3] = [P3_ratio] * T3
								states0[starting + T1 + T2:starting + T1 + T2 + T3] = [3] * T3
								for n in range(starting, ending):
									powratio0[n] = powratio0[n] + 0.05 * P_noise[n - starting]
					timesec, powratio, states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states                    
					Appliance['processed'] = True
					continue
				# ############################## TV ############################
				if name == 'LCDTV' or name == 'CRTTV':
					ending = 0
					for time_step in range(defined_number):
						if Occupancy_pro[time_step] == 1.0 and timesec[60 * time_step] == 0:
							new_index=time_step%(24*60)                            
							threshold = prob_curve[new_index] * n_cycles * k * c
							if random.random() * 100 < threshold:
								duration = min_duration + math.floor((max_duration - min_duration) * random.random())
								starting = int(60 * time_step + math.floor(random.random() * 60))
								ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
								timesec0[starting:ending] = [1.0] * (ending - starting)
								powratio0[starting:ending] = [1.0] * (ending - starting)
								states0[starting:ending] = [1] * (ending - starting)
								# Modify the powratio array by adding noise and falling spike
								n_fallspike = 1 + int(round(2 * random.random()))
								for n in range(n_fallspike):
									t_fall = int(round(starting + (ending - starting) * random.random()))
									powratio0[t_fall:t_fall + 2] = [random.normalvariate(0.2, 0.02)] * 2
								for n in range(starting, ending):
									powratio0[n] = powratio0[n] + 0.02 * P_noise[n - starting]
					timesec, powratio, states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states 
					Appliance['processed'] = True
					continue
				# ############################## Dishwasher ############################
				if name == 'DSW':
					ending = 0
					for time_step in range(defined_number):
						if Occupancy_pro[time_step] == 1.0 and timesec[60 * time_step] == 0:
							new_index=time_step%(24*60) 
							threshold = prob_curve[new_index] * n_cycles * k * c
							if random.random() * 100 < threshold:
								duration = min_duration + math.floor((max_duration - min_duration) * random.random())
								starting = int(60 * time_step + math.floor(random.random() * 60))
								# Modify the powratio array by adding spike and multi-level
								Tfull = int(math.floor(2 * 60 * duration / 3))
								st_dyn = starting
								for n_subcycle in range(2):
									T = int(math.floor(random.normalvariate(240,5)))
									T1 = int(math.floor(random.normalvariate(150,5)))
									T2 = int(math.floor(random.normalvariate(40,2)))
									powratio0[st_dyn:st_dyn + T1] = [random.normalvariate(0.3, 0.003)] * T1
									powratio0[st_dyn + T1:st_dyn + T1 + T2] = [random.normalvariate(0.2, 0.002)] * T2
									powratio0[st_dyn + T1 + T2:st_dyn + T] = [0.01] * (T - T1 - T2)
									states0[st_dyn:st_dyn + T1] = [1] * T1
									states0[st_dyn + T1:st_dyn + T1 + T2] = [2] * T2
									states0[st_dyn + T1 + T2:st_dyn + T] = [3] * (T - T1 - T2)
									Dyn = Dynamic_check(Appliance, name)
									for n_dyn in range(3):
										# When the fidge start each time, insert the 3 sec. dynamic process
										powratio0[st_dyn + n_dyn] = Dyn['Pstrt'][n_dyn]
									st_dyn = st_dyn + T
								powratio0[st_dyn:st_dyn + Tfull] = [1.0] * Tfull
								states0[st_dyn:st_dyn + Tfull] = [4] * Tfull
								st_dyn = st_dyn + Tfull + 40
								for n_subcycle in range(2):
									T = int(math.floor(random.normalvariate(240,5)))
									T1 = int(math.floor(random.normalvariate(150,5)))
									T2 = int(math.floor(random.normalvariate(40,2)))
									powratio0[st_dyn:st_dyn + T1] = [random.normalvariate(0.3, 0.003)] * T1
									powratio0[st_dyn + T1:st_dyn + T1 + T2] = [random.normalvariate(0.2, 0.002)] * T2
									powratio0[st_dyn + T1 + T2:st_dyn + T] = [0.01] * (T - T1 - T2)
									states0[st_dyn:st_dyn + T1] = [1] * T1
									states0[st_dyn + T1:st_dyn + T1 + T2] = [2] * T2
									states0[st_dyn + T1 + T2:st_dyn + T] = [3] * (T - T1 - T2)
									Dyn = Dynamic_check(Appliance, name)
									for n_dyn in range(3):
										# When the fidge start each time, insert the 3 sec. dynamic process
										powratio0[st_dyn + n_dyn] = Dyn['Pstrt'][n_dyn]
									st_dyn = st_dyn + T
								ending = st_dyn
								timesec0[starting:ending] = [1] * (ending - starting)
					timesec, powratio, states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states 
					Appliance['processed'] = True
					continue
				# ############################## PV panel ############################
				if name == 'PV':
					if Appliance['processed'] == True:
						continue
					Pmax = Appliance['P']
					# Calulate the power generation based on the random irradiance
					Power_PV_min = [Pmax * A for A in Irr_pu]
					for n in range(defined_number):
						Power_PV_sec[60 * n:60 * (n + 1)] = [Power_PV_min[n]] * 60
					Appliance['timeofuse'] = [1.0] * n_simulations
					Appliance['powratio'] = Power_PV_sec
					Appliance['processed'] = True
					continue
				# ############################## PHEV panel ############################
				if name == 'PHEV':
					if Appliance['processed'] == True:
						continue
					# Calulate the PHEV usage profile
					time_PHEV_min = PHEVsimulation(House, Appliance)
					for n in range(defined_number):
						time_PHEV_sec[60 * n:60 * (n + 1)] = [time_PHEV_min[n]] * 60
					Appliance['timeofuse'] = time_PHEV_sec
					Appliance['powratio'] = time_PHEV_sec
					Appliance['processed'] = True
					continue
				# ##################### Other usual appliances ########################
				ending = 0
				for time_step in range(defined_number):
					if Occupancy_pro[time_step] == 1.0 and timesec[60 * time_step] == 0:
						new_index=time_step%(24*60) 
						threshold = prob_curve[new_index] * n_cycles * k * c
						if random.random() * 100 < threshold:
							duration = min_duration + math.floor((max_duration - min_duration) * random.random())
							starting = int(60 * time_step + math.floor(random.random() * 60))
							ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
							timesec0[starting:ending] = [1] * (ending - starting)
				powratio0 = timesec0
				states0 = timesec0
				timesec, powratio, states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
				Appliance['timeofuse'] = timesec
				Appliance['powratio'] = powratio
				Appliance['states'] = states 
				Appliance['processed'] = True
			else:
				# SwitchON events are assigned determininstically
				ending = 1
				prob_curve = Appliance['prob_curve']               
				# ########################### PC and LCD ############################
				if name == 'PC' or name == 'LCD':
					if Appliance['processed'] == True:
						continue
					if len(ind_pc) != len(ind_lcd) and enter_pc + 1 > min([len(ind_pc), len(ind_lcd)]):
						ending = 0
						for time_step in range(defined_number):
							new_index=time_step%(24*60) 
							if prob_curve[new_index]:
								duration = min_duration + math.floor((max_duration - min_duration) * random.random())
								starting = int(60 * time_step + math.floor(random.random() * 60))
								ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
								timesec0[starting:ending] = [1.0] * (ending - starting)
								powratio0[starting:ending] = [1.0] * (ending - starting)
								states0[starting:ending] = [1] * (ending - starting)
								# Add multi-stage
								T1 = int(round(60 * duration * (0.2 + 0.4 * random.random())))
								T2 = int(round(60 * duration * (0.2 + 0.5 * random.random())))
								P2_ratio = random.normalvariate(1.3, 0.03)
								powratio0[starting + T1:ending + T1 + T2] = [P2_ratio] * (ending + T2 - starting)
								powratio0[starting + T1:ending + T1 + T2] = [P2_ratio] * (ending + T2 - starting)
								states0[starting + T1:ending + T1 + T2] = [2] * (ending + T2 - starting)
								states0[starting + T1:ending + T1 + T2] = [2] * (ending + T2 - starting)
								for i in range(starting,ending):
									powratio0[i] = powratio0[i] + 0.04 * P_noise[i - starting]
						timesec, powratio, states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
						Appliance['timeofuse'] = timesec
						Appliance['powratio'] = powratio
						Appliance['states'] = states 
						Appliance['pfratio'] = pfratio
						Appliance['processed'] = True
						enter_pc = enter_pc + 1
					else:
						ending = 0
						for time_step in range(defined_number):
							if prob_curve[time_step]:
								# PC
								duration = min_duration + math.floor((max_duration - min_duration) * random.random())
								starting = int(60 * time_step + math.floor(random.random() * 60))
								ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
								timesec0[starting:ending] = [1.0] * (ending - starting)
								powratio0[starting:ending] = [1.0] * (ending - starting)
								states0[starting:ending] = [1] * (ending - starting)
								# Add multi-stage
								T1 = int(round(duration * (0.2 + 0.4 * random.random())))
								T2 = int(round(duration * (0.2 + 0.5 * random.random())))
								P2_ratio = random.normalvariate(1.3, 0.03)
								powratio0[starting + T1:ending + T1 + T2] = [P2_ratio] * (ending + T2 - starting)
								powratio0[starting + T1:ending + T1 + T2] = [P2_ratio] * (ending + T2 - starting)
								states0[starting + T1:ending + T1 + T2] = [2] * (ending + T2 - starting)
								states0[starting + T1:ending + T1 + T2] = [2] * (ending + T2 - starting)
								for i in range(starting,ending):
									powratio0[i] = powratio0[i] + 0.08 * P_noise[i - starting]
						timesec, powratio, states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
						House['appliances'][ind_pc[enter_pc] - 1]['timeofuse'] = timesec
						House['appliances'][ind_pc[enter_pc] - 1]['powratio'] = powratio
						House['appliances'][ind_pc[enter_pc] - 1]['states'] = states
						House['appliances'][ind_pc[enter_pc] - 1]['processed'] = True
						House['appliances'][ind_lcd[enter_pc] - 1]['timeofuse'] = timesec
						House['appliances'][ind_lcd[enter_pc] - 1]['powratio'] = powratio
						House['appliances'][ind_lcd[enter_pc] - 1]['states'] = states
						House['appliances'][ind_lcd[enter_pc] - 1]['processed'] = True
						enter_pc = enter_pc + 1
					continue
				# ################ Washer, dryer and ASD dryer ##########################
				if name == 'WSH' or name == 'DRY' or name == 'ASDDRY':
					if Appliance['processed'] == True:
						continue
					min_duration_wsh = app.appliance_char[House['appliances'][ind_wsh[enter_wshdry] - 1]['app_name']]['min_cycleduration']
					max_duration_wsh = app.appliance_char[House['appliances'][ind_wsh[enter_wshdry] - 1]['app_name']]['max_cycleduration']
					min_duration_dry = app.appliance_char[House['appliances'][ind_dry[enter_wshdry] - 1]['app_name']]['min_cycleduration']
					max_duration_dry = app.appliance_char[House['appliances'][ind_dry[enter_wshdry] - 1]['app_name']]['max_cycleduration']
					ending_wsh = 0
					ending_dry = 0
					for time_step in range(defined_number):
						new_index=time_step%(24*60)                        
						if prob_curve[new_index]:
							# Washing machine
							duration = min_duration + math.ceil((max_duration_wsh - min_duration_wsh) * random.random())
							starting_wsh = int(60 * time_step + math.floor(random.random() * 60))
							starting_wsh0 = starting_wsh
							ending_wsh = int(60 * (time_step + duration) + math.floor(random.random() * 60))
							timesec_wsh0[starting_wsh:ending_wsh] = [1.0] * (ending_wsh - starting_wsh)
							# Add pulse
							while True:
								T = int(round(random.normalvariate(30,1)))
								Ton = int(T * round(0.5 + 0.1 * random.random())) # The duration of motor-on
								powratio_wsh0[starting_wsh:starting_wsh + Ton] = [1.0] * Ton
								powratio_wsh0[starting_wsh + Ton:starting_wsh + Ton + T] = [0.01] * Ton
								states_wsh0[starting_wsh:starting_wsh + Ton] = [1] * Ton
								states_wsh0[starting_wsh + Ton:starting_wsh + Ton + T] = [2] * Ton
								starting_wsh = starting_wsh + T
								if starting_wsh > ending_wsh:
									break
							# Add noise
							for i in range(starting_wsh0, ending_wsh):
								powratio_wsh0[i] = powratio_wsh0[i] + 0.06 * P_noise[i - starting_wsh0]
							# dryer
							duration = int(math.floor(min_duration_dry + (max_duration_dry - min_duration_dry) * random.random()))
							starting_dry = int(60 * (time_step + 40) + math.floor(60 * random.random()))
							starting_dry0 = starting_dry
							ending_dry = int(60 * (time_step + duration + 40) + math.floor(random.random() * 60))
							timesec0[starting_dry:ending_dry] = [1.0] * (ending_dry - starting_dry)
							for n in range(starting_dry, ending_dry):
								powratio0[n] = powratio0[n] + (-0.05 + 0.06 * random.random())
								states0[n] = states0[n] + (1)
							while True:
								T = int(round(random.normalvariate(240,20)))
								Ton = int(T * round(0.5 + 0.1 * random.random())) # The duration of motor-on
								for n_decr in range(10):
									powratio0[starting_dry + n_decr] = 1.07 - 0.007 * (n_decr - 1)
									states0[starting_dry + n_decr] =2
								powratio0[starting_dry + 10:starting_dry + Ton] = [1.0] * (Ton - 10)
								powratio0[starting_dry + Ton:starting_dry + T] = [0.001] * (T - Ton)
								states0[starting_dry + 10:starting_dry + Ton] = [3] * (Ton - 10)
								states0[starting_dry + Ton:starting_dry + T] = [4] * (T - Ton)
								starting_dry = starting_dry + T
								if starting_dry > ending_dry:
										break
							for i in range(starting_dry0, ending_dry):
								powratio0[i] = powratio0[i] + 0.005 * P_noise[i - starting_dry0]
					# Check if the operation time goes to the next day
					timesec_wsh, powratio_wsh,states_wsh = TimeModify2(timesec_wsh0, powratio_wsh0,states_wsh0,ending_wsh)
					# Check if the operation time goes to the next day
					timesec, powratio, states = TimeModify2(timesec0, powratio0, states0, ending_dry)
					# Add duty cycle to the DRY, 1 min per small cycle
					House['appliances'][ind_wsh[enter_wshdry] - 1]['timeofuse'] = timesec_wsh
					House['appliances'][ind_wsh[enter_wshdry] - 1]['powratio'] = powratio_wsh
					House['appliances'][ind_wsh[enter_wshdry] - 1]['states'] = states_wsh
					House['appliances'][ind_wsh[enter_wshdry] - 1]['processed'] = True
					House['appliances'][ind_dry[enter_wshdry] - 1]['timeofuse'] = timesec
					House['appliances'][ind_dry[enter_wshdry] - 1]['powratio'] = powratio
					House['appliances'][ind_dry[enter_wshdry] - 1]['states'] = states
					House['appliances'][ind_dry[enter_wshdry] - 1]['processed'] = True
					enter_wshdry = enter_wshdry + 1
					continue
				# ########################## Stove ################################
				if name == 'STO':
					ending = 0
					for time_step in range(defined_number):
						new_index=time_step%(24*60)                         
						if prob_curve[new_index]:
							duration = min_duration + math.ceil((max_duration - min_duration) * random.random())
							starting = int(60 * time_step + math.floor(random.random() * 60))
							ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
							timesec0[starting:ending] = [1.0] * (ending - starting)
							Ton_ratio = 0.15 + 0.25 * random.random()		# The duration of motor-on
							while True:
								T = int(round(random.normalvariate(30,2)))
								Ton = int(round(Ton_ratio * T))
								powratio0[starting:starting + Ton] = [1.0] * Ton
								powratio0[starting + Ton:starting + T] = [0.001] * (T - Ton)
								states0[starting:starting + Ton] = [1] * Ton
								states0[starting + Ton:starting + T] = [2] * (T - Ton)
								starting = starting + T
								if starting > ending:
									break
					# Check if the operation time goes to the next
					timesec, powratio,states = TimeModify2(timesec0, powratio0, states0, ending)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states
					Appliance['processed'] = True
					continue
				# ########################## Coffee maker, heater ##########################
				if name == 'COF' or name == 'HEA':
					ending = 0
					for time_step in range(defined_number):
						new_index=time_step%(24*60)                         
						if prob_curve[new_index]:
							duration = min_duration + math.ceil((max_duration - min_duration) * random.random())
							starting = int(60 * time_step + math.floor(random.random() * 60))
							ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
							timesec0[starting:ending] = [1.0] * (ending - starting)
							if name == 'COF':
								T = int(round(random.normalvariate(30,2)))
								Ton = int(round(T * (0.3 + 0.4 * random.random())))
								while True:
									powratio0[starting:starting + Ton] = [1.0] * Ton
									powratio0[starting + Ton:starting + T] = [0.001] * (T - Ton)
									states0[starting:starting + Ton] = [1] * Ton
									states0[starting + Ton:starting + T] = [2] * (T - Ton)
									starting = starting + T
									if starting > ending:
										break
							T = int(round(random.normalvariate(40,2)))
							Ton = int(round(T * (0.1 + 0.4 * random.random())))
							while True:
								powratio0[starting:starting + Ton] = [1.0] * Ton
								powratio0[starting + Ton:starting + T] = [0.01] * (T - Ton)
								states0[starting:starting + Ton] = [1] * Ton
								states0[starting + Ton:starting + T] = [2] * (T - Ton)
								starting = starting + T
								if starting > ending:
									break
					timesec, powratio,states = TimeModify2(timesec0, powratio0,states0, ending)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states
					Appliance['processed'] = True
					continue
				# ########################## Vacuum Cleaner ##########################
				if name == 'VAC':
					for time_step in range(defined_number):
						new_index=time_step%(24*60)                         
						if prob_curve[new_index]:
							duration = min_duration + math.ceil((max_duration - min_duration) * random.random())
							starting = int(60 * time_step + math.floor(random.random() * 60))
							ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
							timesec0[starting:ending] = [1.0] * (ending - starting)
							powratio0 = timesec0
							powratio0[starting] = 1.5 - 0.25 * 2
							states0 = timesec0
							states0[starting] = 2
					timesec, powratio,states = TimeModify2(timesec0, powratio0,states0, ending)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states
					Appliance['processed'] = True
					continue
				# ########################## Laptop ##########################
				if name == 'LAP':
					ending = 0
					for time_step in range(defined_number):
						new_index=time_step%(24*60)                         
						if prob_curve[new_index]:
							duration = min_duration + math.floor((max_duration - min_duration) * random.random())
							starting = int(60 * time_step + math.floor(random.random() * 60))
							ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
							timesec0[starting:ending] = [1.0] * (ending - starting)
							# Add multi-stage
							T1 = max([int(round(random.normalvariate(30,5))), 10])
							P1_ratio = max([random.normalvariate(1.4,0.05), 1.2])
							powratio0[starting:starting + T1] = [P1_ratio] * T1
							powratio0[starting + T1:ending] = [1.0] * (ending - starting - T1)
							states0[starting:starting + T1] = [1] * T1
							states0[starting + T1:ending] = [2] * (ending - starting - T1)
							T2 = int(round(60 * duration * (0.1 + 0.4 * random.random())))
							T3 = int(round(60 * duration * (0.2 + 0.5 * random.random())))
							P3_ratio = random.normalvariate(1.3,0.05)
							powratio0[starting + T1 + T2:starting + T1 + T2 + T3] = [P3_ratio] * T3
							states0[starting + T1 + T2:starting + T1 + T2 + T3] = [3] * T3
							for n in range(starting, ending):
								powratio0[n] = powratio0[n] + 0.05 * P_noise[n - starting]
					timesec, powratio,states = TimeModify2(timesec0, powratio0,states0, ending)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states
					Appliance['processed'] = True
					continue
				# ############################## TV ############################
				if name == 'LCDTV' or name == 'CRTTV':
					ending = 0
					for time_step in range(defined_number):
						new_index=time_step%(24*60)                         
						if prob_curve[new_index]:
							duration = min_duration + math.floor((max_duration - min_duration) * random.random())
							starting = int(60 * time_step + math.floor(random.random() * 60))
							ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
							timesec0[starting:ending] = [1.0] * (ending - starting)
							powratio0 = timesec0
							states0 = timesec0
							# Modify the powratio array by adding noise and falling spike
							n_fallspike = 1 + int(round(2 * random.random()))
							for n in range(n_fallspike):
								t_fall = int(round(starting + (ending - starting) * random.random()))
								powratio0[t_fall:t_fall + 2] = [random.normalvariate(0.2, 0.02)] * 2
								sttates0[t_fall:t_fall + 2] = [2] * 2
							for n in range(starting, ending):
								powratio0[n] = powratio0[n] + 0.02 * P_noise[n - starting]
					timesec, powratio,states = TimeModify2(timesec0, powratio0,states0, ending)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states
					Appliance['processed'] = True
					continue
				# ############################## Dishwasher ############################
				if name == 'DSW':
					ending = 0
					for time_step in range(defined_number):
						new_index=time_step%(24*60)                         
						if prob_curve[new_index]:
							duration = min_duration + math.floor((max_duration - min_duration) * random.random())
							starting = int(60 * time_step + math.floor(random.random() * 60))
							# Modify the powratio array by adding spike and multi-level
							Tfull = int(math.floor(2 * 60 * duration / 3))
							st_dyn = starting
							for n_subcycle in range(2):
								T = int(math.floor(random.normalvariate(240,5)))
								T1 = int(math.floor(random.normalvariate(150,5)))
								T2 = int(math.floor(random.normalvariate(40,2)))
								powratio0[st_dyn:st_dyn + T1] = [random.normalvariate(0.3, 0.003)] * T1
								powratio0[st_dyn + T1:st_dyn + T1 + T2] = [random.normalvariate(0.2, 0.002)] * T2
								powratio0[st_dyn + T1 + T2:st_dyn + T] = [0.01] * (T - T1 - T2)
								states0[st_dyn:st_dyn + T1] = [1] * T1
								states0[st_dyn + T1:st_dyn + T1 + T2] = [2] * T2
								states0[st_dyn + T1 + T2:st_dyn + T] = [3] * (T - T1 - T2)
								Dyn = Dynamic_check(Appliance, name)
								for n_dyn in range(3):
									# When the fidge start each time, insert the 3 sec. dynamic process
									powratio0[st_dyn + n_dyn] = Dyn['Pstrt'][n_dyn]
								st_dyn = st_dyn + T
							powratio0[st_dyn:st_dyn + Tfull] = [1.0] * Tfull
							states0[st_dyn:st_dyn + Tfull] = [4] * Tfull
							st_dyn = st_dyn + Tfull + 40
							for n_subcycle in range(2):
								T = int(math.floor(random.normalvariate(240,5)))
								T1 = int(math.floor(random.normalvariate(150,5)))
								T2 = int(math.floor(random.normalvariate(40,2)))
								powratio0[st_dyn:st_dyn + T1] = [random.normalvariate(0.3, 0.003)] * T1
								powratio0[st_dyn + T1:st_dyn + T1 + T2] = [random.normalvariate(0.2, 0.002)] * T2
								powratio0[st_dyn + T1 + T2:st_dyn + T] = [0.01] * (T - T1 - T2)
								states0[st_dyn:st_dyn + T1] = [1] * T1
								states0[st_dyn + T1:st_dyn + T1 + T2] = [2] * T2
								states0[st_dyn + T1 + T2:st_dyn + T] = [3] * (T - T1 - T2)
								Dyn = Dynamic_check(Appliance, name)
								for n_dyn in range(3):
									# When the fidge start each time, insert the 3 sec. dynamic process
									powratio0[st_dyn + n_dyn] = Dyn['Pstrt'][n_dyn]
								st_dyn = st_dyn + T
							ending = st_dyn
							timesec0[starting:ending] = [1] * (ending - starting)
					timesec, powratio, states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states
					Appliance['processed'] = True
					continue
				# ############################# Treadmill #############################
				if name == 'TRD':
					ending = 0
					for time_step in range(defined_number):
						new_index=time_step%(24*60)                         
						if prob_curve[new_index]:
							duration = min_duration + math.floor((max_duration - min_duration) * random.random())
							starting = int(60 * time_step + math.floor(random.random() * 60))
							ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
							timesec0[starting:ending] = [1] * (ending - starting)
							powratio0[starting:ending] = [1] * (ending - starting)
							states0[starting:ending] = [1] * (ending - starting)
							for n in range(starting, ending):
								powratio0[n] = powratio0[n] + 0.04 * P_noise[n - starting]
					timesec, powratio,states = TimeModify2(timesec0, powratio0, states0,ending,n_simulations)
					Appliance['timeofuse'] = timesec
					Appliance['powratio'] = powratio
					Appliance['states'] = states
					Appliance['processed'] = True
					continue
				# ##################### Other usual appliances ########################
				ending = 0
				for time_step in range(defined_number):
					new_index=time_step%(24*60)                     
					if prob_curve[new_index]:
						duration = min_duration + math.floor((max_duration - min_duration) * random.random())
						starting = int(60 * time_step + math.floor(random.random() * 60))
						ending = int(60 * (time_step + duration) + math.floor(random.random() * 60))
						timesec0[starting:ending] = [1] * (ending - starting)
				powratio0 = timesec0
				states0 = timesec0
				timesec, powratio,states = TimeModify2(timesec0, powratio0,states0, ending)
				Appliance['timeofuse'] = timesec
				Appliance['powratio'] = powratio
				Appliance['states'] = states
				Appliance['processed'] = True
	return app

def PVsimulation(day_type, Irr_full,defined_number):
	import PV_database
	import math
	import random
	if day_type == 1:
		Cloud_block_level = PV_database.Cloud_block['Clr']
	elif day_type == 2:
		Cloud_block_level = PV_database.Cloud_block['Most_clear']
	elif day_type == 3:
		Cloud_block_level = PV_database.Cloud_block['Mostly_cloudy']
	elif day_type == 4:
		Cloud_block_level = PV_database.Cloud_block['Overcast']
	CCL_corr = [0] * defined_number
	Irr_pu = [0] * defined_number
	endtime = 0
	count = 1
	for time_step in range(defined_number):
		# If this cloud cycle has not come to an end, we don't have to generate a new CCL (Cloud Coverage Level)
		if time_step + 1 > int(math.ceil(endtime)):
			count = count + 1
			# Determine cloud type
			if day_type == 4:
				threshold2 = Cloud_block_level[1][2]
				aux = random.random() * 100
				if aux < threshold2:
					cover_min = Cloud_block_level[1][0]
					cover_max = Cloud_block_level[1][1]
				else:
					cover_min = Cloud_block_level[0][0]
					cover_max = Cloud_block_level[0][1]
			else:
				threshold2 = Cloud_block_level[1][2] + Cloud_block_level[2][2]
				threshold3 = Cloud_block_level[2][2]
				aux = random.random() * 100
				if aux < threshold3:
					cover_min = Cloud_block_level[2][0]
					cover_max = Cloud_block_level[2][1]
				elif aux < threshold2:
					cover_min = Cloud_block_level[1][0]
					cover_max = Cloud_block_level[1][1]
				else:
					cover_min = Cloud_block_level[0][0]
					cover_max = Cloud_block_level[0][1]
			# Determine Cloud Coverage Level(CCL)
			CCL = (cover_min + (cover_max - cover_min) * random.random())/100
			# Determine Cloud Cycle (endtime)
			duration = 24 * random.random()
			duration = duration * 60
			endtime = time_step + duration
			if endtime > defined_number:
				endtime = defined_number
		# CCL correction
		CCL_corr[time_step] = random.normalvariate(CCL, 0.04)
		
		mod_time_step=(time_step)% (int(3600 * 24 / 60))
		
		Irr_pu[time_step] = PV_database.Max_irradiance[mod_time_step] * (1 - CCL_corr[time_step])/Irr_full
		# Prevent the PV from accepting irradiance larger than full irradiance
		if Irr_pu[time_step] > 1:
			Irr_pu[time_step] = 1
		# Prevent the PV from generating negative random power
		if Irr_pu[time_step] < 0:
			Irr_pu[time_step] = 0
	return Irr_pu

def PHEVsimulation(House, Appliance):
	import PHEV_database	# Load the PHEV related data
	import random
	import math
	time_PHEV = [0] * (60 * 24)
	# Arrange appliance working time
	char_strategy = PHEV_database.strategy
	travel = PHEV_database.travel['default_pro']
	Power = Appliance['P']
	if House['PHEV_strategy'] == 'uncontrolled':
		charging_strategy = 0
	elif House['PHEV_strategy'] == 'controlled':
		charging_strategy = 1
	elif House['PHEV_strategy'] == 'smart':
		charging_strategy = 2
	# ################## PHEV charging start time ###################
	PHEV_cdf = char_strategy[charging_strategy]
	for i in range(1,len(PHEV_cdf)):
		PHEV_cdf[i] = PHEV_cdf[i-1] + PHEV_cdf[i]
	start_cdf = random.random() * 100
	# We assume that the PHEV is charged only once a day
	for n in range(defined_number):
		if start_cdf < PHEV_cdf[n]:
			chraging_start = int(n)
			break
	if n == 1439:
		chraging_start = 1200
	# ################# PHEV charging duration #####################
	Ecar_range = 40.0		# All electric range (maximal miles the PHEV can run with battery only)
	Ecar_capacity = PHEV_database.battery['capacity']	# Battery capacity is 16 kW, with efficiency 0.88
	DOD = PHEV_database.battery['DOD']
	yita = PHEV_database.battery['efficiency']
	for i in range(1,len(travel)):
		travel[i] = travel[i-1] + travel[i]
	Travel_cdf = travel
	start_travel_cdf = random.random() * 100
	for n in range(1001):
		if start_travel_cdf < Travel_cdf[n]:
			travel_mile = (n / 10.0)
			break
	if travel_mile > 40:
		SOC_battery = 0
	else:
		SOC_battery = (Ecar_range - travel_mile) / Ecar_range
	charging_duration = int(60 * math.ceil(1000 * Ecar_capacity * (1 - SOC_battery) * DOD / (Power * yita))) # 60 * (1000 * kWh/W), the unit of charging duration is min.
	time_PHEV[chraging_start:chraging_start + charging_duration] = [1.0] * charging_duration
	len_time_PHEV = len(time_PHEV)
	if len_time_PHEV > defined_number:
		time_PHEV[:len_time_PHEV - defined_number] = time_PHEV[1441:len_time_PHEV]
		time_PHEV = time_PHEV[:defined_number]
	return time_PHEV

def TimeModify2(timesec0, powratio0,states0, ending, n_simulations):
	timesec = [0] * n_simulations
	powratio = [0] * n_simulations
	states=[0] * n_simulations
	if ending > n_simulations:
		timesec[:] = timesec0[:n_simulations]
		timesec[:ending - n_simulations] = timesec0[n_simulations:ending]
		powratio[:] = powratio0[:n_simulations]
		powratio[:ending - n_simulations] = powratio0[n_simulations:ending]
		states[:] = states0[:n_simulations]
		states[:ending - n_simulations] = states0[n_simulations:ending]
	else:
		timesec = timesec0[0:n_simulations]
		powratio = powratio0[0:n_simulations]
		states = states0[0:n_simulations]
	return timesec, powratio, states

def TimeModify(timesec0, powratio0, ending, n_simulations):
	timesec = [0] * n_simulations
	powratio = [0] * n_simulations
	if ending > n_simulations:
		timesec[:] = timesec0[:n_simulations]
		timesec[:ending - n_simulations] = timesec0[n_simulations:ending]
		powratio[:] = powratio0[:n_simulations]
		powratio[:ending - n_simulations] = powratio0[n_simulations:ending]
	else:
		timesec = timesec0[0:n_simulations]
		powratio = powratio0[0:n_simulations]
	return timesec, powratio

def Dynamic_check(Appliance, name):
	import math
	if name == 'FRZR' or name == 'RFR' or name == 'ASDFR':
		Papp = Appliance['app_qtt'] * Appliance['P']
		PFapp = Appliance['pf']
		Pmotor_in = 0.8 * Papp
		Pm, PFm = Motorstart_cmprs(Pmotor_in, PFapp)
		Pstrt = [A / Papp + 0.2 for A in Pm]
		Dyn = {}
		Dyn['Pstrt'] = Pstrt
		Dyn['PFstrt'] = PFm
	elif name == 'FUR':
		Papp = Appliance['app_qtt'] * Appliance['P']
		PFapp = Appliance['pf']
		Pmotor_in = 0.2 * Papp
		Pm, PFm = Motorstart_fan(Pmotor_in, PFapp)
		Pstrt = [A / Papp + 0.8 for A in Pm]
		Qstrt = [0] * 3
		PFstrt = [0] * 3
		for j in range(3):
			Qstrt[j] = 0.8 * Papp * math.sqrt(1 - PFapp**2)/PFapp + Pm[j] * math.sqrt(1 - PFm[j]**2)/PFm[j]
			PFstrt[j] = Pstrt[j] / (math.sqrt(Pstrt[j]**2 + Qstrt[j]**2))
		Dyn = {}
		Dyn['Pstrt'] = Pstrt
		Dyn['PFstrt'] = PFstrt
	elif name == 'DSW':
		Papp = Appliance['app_qtt'] * Appliance['P']
		PFapp = Appliance['pf']
		Pmotor_in = 0.1 * Papp
		Pm, PFm = Motorstart_cmprs(Pmotor_in, PFapp)
		Pstrt = [A / Papp + 0.2 for A in Pm]
		Dyn = {}
		Dyn['Pstrt'] = Pstrt
		Dyn['PFstrt'] = PFm
	return Dyn

def Motorstart_cmprs(Pmotor_in, PFapp):
	# Preprocess the data
	import math
	import random
	delay = random.randint(1, 29)
	Vs = 120			# Voltage
	Pmotor_mec = Pmotor_in * (0.5 + 0.1 * random.random()) # Mechanical power = total_power * motor_ratio * efficiency (0.5 ~ 0.7)
	R1m_base = 6.0
	X1m_base = 8.0
	R2m_base = 15.0/2.0
	X2m_base = 10.0/2.0
	R1a_base = 6.0
	X1a_base = 8.0
	Xem_base = 150.0/2.0
	Zs = 0.02  # Numeros y mas numeros magicos!!
	
	a = 0.65	# The ratio between auxiliary turns and main turns
	rotating_s = 1800	# rpm
	omega_s = 2 * 3.1416 * rotating_s / 60 # rad/sec.
	slip_steady = 0.04 + 0.01 * random.random()
	Tload_base = Pmotor_mec / (omega_s * (1 - slip_steady)) # N*m
	J_base = 0.01	# kg * m**2
	N_per_sec = 30
	dT = 1.0 / N_per_sec	# sec.
	len = 150	# number of snapshots
	
	# Search for the optimal running capacitor
	Xcr = [0] * 25
	PF_scan = [0] * 25
	for j in range(25):
		Cr = 25 + j
		Xcr[j] = 1e6 / (377 * Cr)
		Zf_scan = 1j * Xem_base * (R2m_base / slip_steady + 1j * X2m_base) / (R2m_base / slip_steady + 1j * (Xem_base + X2m_base))
		Zb_scan = 1j * Xem_base * (R2m_base / (2 - slip_steady) + 1j * X2m_base) / (R2m_base / (2 - slip_steady) + 1j * (Xem_base + X2m_base))
		Zm_scan = R1m_base + 1j * X1m_base + Zf_scan + Zb_scan
		Za_scan = R1a_base + 1j * X1a_base + a**2 * (Zf_scan + Zb_scan) - 1j * Xcr[j]
		Ztotal_scan = Zm_scan * Za_scan / (Zm_scan + Za_scan)
		PF_scan[j] = math.cos(math.atan2(Ztotal_scan.imag, Ztotal_scan.real))
		if PF_scan[j] > PFapp - 0.02 and PF_scan[j] < PFapp + 0.02:
			break
	Xcr_base = Xcr[j]

	# Transform the bse motor parameters into realistic motor parameters
	P_pu = Pmotor_in/100.0
	R1m = R1m_base/P_pu
	X1m = X1m_base/P_pu
	R2m = R2m_base/P_pu
	X2m = X2m_base/P_pu
	R1a = R1a_base/P_pu
	X1a = X1a_base/P_pu
	Xem = Xem_base/P_pu
	Xcr = Xcr_base/P_pu
	J = J_base * P_pu
	Tload = Tload_base * P_pu
	
	# Optimize the starting capacitor to obtain the maximum starting torque per Ampere
	Zf = 1j * Xem * (R2m + 1j * X2m) / (R2m + 1j * (Xem + X2m))
	Zb = 1j * Xem * (R2m + 1j * X2m) / (R2m + 1j * (Xem + X2m));
	Zm = R1m + 1j * X1m + Zf + Zb
	Za_noc = R1a + 1j * X1a + a**2 * (Zf + Zb)
	Ra = Za_noc.real 
	Xa = Za_noc.imag 
	Rm = Zm.real
	Xm = Zm.imag
	Xc = Xa + (-Xm * Ra + math.sqrt(Rm**2 + Xm**2) * math.sqrt(Ra * (Ra + Rm)))/Rm;

	# Calculate the dynamic parameters through iteration
	Im			= [0] * len
	Ia			= [0] * len
	I			= [0] * len
	Vrec		= [0] * len
	slip		= [0] * len
	torque_mec	= [0] * len
	omega_m		= [0] * len
	Ztotal		= [0] * len
	VIAng		= [0] * len
	PF			= [0] * len
	Pin			= [0] * len
	omega_m[0]	= 2 * 3.1416 * 0 / 60 # Ojo!! numero magico
	for k in range(len - 1):
		if omega_m[k] < 0.75 * omega_s:
			# Auxiliary winding is connected
			slip[k] = (omega_s - omega_m[k]) / omega_s
			Zf = 1j * Xem * (R2m / slip[k] + 1j * X2m) / (R2m / slip[k] + 1j * (Xem + X2m))
			Zb = 1j * Xem * (R2m / (2 - slip[k]) + 1j * X2m) / (R2m / (2 - slip[k]) + 1j * (Xem + X2m))
			Rf = Zf.real 
			Rb = Zb.real
			
			Zm = R1m + 1j * X1m + Zf + Zb
			Za = R1a + 1j * X1a + a**2 * (Zf + Zb) - 1j * Xc
			Ztotal = Zm * Za / (Zm + Za)
			I[k] = Vs / (Ztotal + Zs)
			Vrec[k] = Vs - I[k] * Zs
			Im[k] = Vrec[k] / Zm
			Ia[k] = Vrec[k] / Za
			
			theta = math.atan2(Ia[k].imag, Ia[k].real) - math.atan2(Im[k].imag, Im[k].real)
			torque_mec[k] = 2 * a * math.sqrt(Im[k].real**2 + Im[k].imag**2) * math.sqrt(Ia[k].real**2 + Ia[k].imag**2) * (Rb + Rf) * math.sin(theta) / omega_s
			omega_m[k + 1] = (dT / J) * (torque_mec[k] - Tload) + omega_m[k]
			
			VIAng[k] = math.atan2(Vrec[k].imag, Vrec[k].real) - math.atan2(I[k].imag, I[k].real)
			PF[k] = math.cos(VIAng[k])
			Pin[k] = math.sqrt(Vrec[k].real**2 + Vrec[k].imag**2) * math.sqrt(I[k].real**2 + I[k].imag**2) * PF[k]
		else:
			# When the rotating speed is large enough, the auxiliary winding is disconnected      
			slip[k] = (omega_s - omega_m[k]) / omega_s
			Zf = 1j * Xem * (R2m / slip[k] + 1j * X2m) / (R2m / slip[k] + 1j * (Xem + X2m))
			Zb = 1j * Xem * (R2m / (2 - slip[k]) + 1j * X2m) / (R2m / (2 - slip[k]) + 1j * (Xem + X2m))
			Rf = Zf.real
			Rb = Zb.real
			
			# Capacitor - start type
			Zm = R1m + 1j * X1m + Zf + Zb
			Za = R1a + 1j * X1a + a**2 * (Zf + Zb) - 1j * Xcr
			Ztotal = Zm * Za / (Zm + Za)
			I[k] = Vs / (Ztotal + Zs)
			Vrec[k] = Vs - I[k] * Zs
			Im[k] = Vrec[k] / Zm
			Ia[k] = Vrec[k] / Za
			
			theta = math.atan2(Ia[k].imag, Ia[k].real) - math.atan2(Im[k].imag, Im[k].real)
			torque_mec[k] = 2 * a * math.sqrt(Im[k].real**2 + Im[k].imag**2) * math.sqrt(Ia[k].real**2 + Ia[k].imag**2) * (Rb + Rf) * math.sin(theta) / omega_s
			omega_m[k + 1] = (dT / J) * (torque_mec[k] - Tload) + omega_m[k]
			
			VIAng[k] = math.atan2(Vrec[k].imag, Vrec[k].real) - math.atan2(I[k].imag, I[k].real)
			PF[k] = math.cos(VIAng[k])
			Pin[k] = math.sqrt(Vrec[k].real**2 + Vrec[k].imag**2) * math.sqrt(I[k].real**2 + I[k].imag**2) * PF[k]
	# Transform the dynamic power parameters into equivalent impedance for load flow calculation
	Pstrt = [0] * 3
	PFstrt = [0] * 3
	for j in range(3):
		index = 30 * j + 2 - delay
		if index < 1:
			Pstrt[j] = 1
			PFstrt[j] = PFapp
		else:
			Pstrt[j] = 0.9 * Pin[index - 1]
			PFstrt[j] = PFapp
	return Pstrt, PFstrt

def Motorstart_fan(Pmotor_in, PFapp):

	# Preprocess the data
	import math
	import random
	delay = random.randint(1, 29)
	Vs = 120			# Voltage
	Pmotor_mec = Pmotor_in * (0.5 + 0.1 * random.random()) # Mechanical power = total_power * motor_ratio * efficiency (0.5 ~ 0.7)
	R1m_base = 6.0
	X1m_base = 8.0
	R2m_base = 15.0/2.0
	X2m_base = 10.0/2.0
	R1a_base = 6.0
	X1a_base = 8.0
	Xem_base = 150.0/2.0
	Zs = 0.02  # Numeros y mas numeros magicos!!
	
	a = 0.65	# The ratio between auxiliary turns and main turns
	rotating_s = 1800	# rpm
	omega_s = 2 * 3.1416 * rotating_s / 60 # rad/sec.
	slip_steady = 0.04 + 0.01 * random.random()
	Tload_base = Pmotor_mec / (omega_s * (1 - slip_steady)) # N*m
	J_base = 0.01	# kg * m**2
	N_per_sec = 30
	dT = 1.0 / N_per_sec	# sec.
	len = 150	# number of snapshots
	
	# Search for the optimal running capacitor
	Xcr = [0] * 25
	PF_scan = [0] * 25
	for j in range(25):
		Cr = 25 + j
		Xcr[j] = 1e6 / (377 * Cr)
		Zf_scan = 1j * Xem_base * (R2m_base / slip_steady + 1j * X2m_base) / (R2m_base / slip_steady + 1j * (Xem_base + X2m_base))
		Zb_scan = 1j * Xem_base * (R2m_base / (2 - slip_steady) + 1j * X2m_base) / (R2m_base / (2 - slip_steady) + 1j * (Xem_base + X2m_base))
		Zm_scan = R1m_base + 1j * X1m_base + Zf_scan + Zb_scan
		Za_scan = R1a_base + 1j * X1a_base + a**2 * (Zf_scan + Zb_scan) - 1j * Xcr[j]
		Ztotal_scan = Zm_scan * Za_scan / (Zm_scan + Za_scan)
		PF_scan[j] = math.cos(math.atan2(Ztotal_scan.imag, Ztotal_scan.real))
		if PF_scan[j] > PFapp - 0.02 and PF_scan[j] < PFapp + 0.02:
			break
	Xcr_base = Xcr[j]

	# Transform the bse motor parameters into realistic motor parameters
	P_pu = Pmotor_in/100.0
	R1m = R1m_base/P_pu
	X1m = X1m_base/P_pu
	R2m = R2m_base/P_pu
	X2m = X2m_base/P_pu
	R1a = R1a_base/P_pu
	X1a = X1a_base/P_pu
	Xem = Xem_base/P_pu
	Xcr = Xcr_base/P_pu
	J = J_base * P_pu
	Tload = Tload_base * P_pu
	
	# Optimize the starting capacitor to obtain the maximum starting torque per Ampere
	Zf = 1j * Xem * (R2m + 1j * X2m) / (R2m + 1j * (Xem + X2m))
	Zb = 1j * Xem * (R2m + 1j * X2m) / (R2m + 1j * (Xem + X2m));
	Zm = R1m + 1j * X1m + Zf + Zb
	Za_noc = R1a + 1j * X1a + a**2 * (Zf + Zb)
	Ra = Za_noc.real 
	Xa = Za_noc.imag 
	Rm = Zm.real
	Xm = Zm.imag
	Xc = Xa + (-Xm * Ra + math.sqrt(Rm**2 + Xm**2) * math.sqrt(Ra * (Ra + Rm)))/Rm;

	# Calculate the dynamic parameters through iteration
	Im			= [0] * len
	Ia			= [0] * len
	I			= [0] * len
	Vrec		= [0] * len
	slip		= [0] * len
	torque_mec	= [0] * len
	omega_m		= [0] * len
	Tload_dy	= [0] * len
	Ztotal		= [0] * len
	VIAng		= [0] * len
	PF			= [0] * len
	Pin			= [0] * len
	omega_m[0]	= 2 * 3.1416 * 0 / 60 # Ojo!! numero magico
	for k in range(len - 1):
		if omega_m[k] < 0.75 * omega_s:
			# Auxiliary winding is connected
			slip[k] = (omega_s - omega_m[k]) / omega_s
			Zf = 1j * Xem * (R2m / slip[k] + 1j * X2m) / (R2m / slip[k] + 1j * (Xem + X2m))
			Zb = 1j * Xem * (R2m / (2 - slip[k]) + 1j * X2m) / (R2m / (2 - slip[k]) + 1j * (Xem + X2m))
			Rf = Zf.real 
			Rb = Zb.real
			
			Zm = R1m + 1j * X1m + Zf + Zb
			Za = R1a + 1j * X1a + a**2 * (Zf + Zb) - 1j * Xc
			Ztotal = Zm * Za / (Zm + Za)
			I[k] = Vs / (Ztotal + Zs)
			Vrec[k] = Vs - I[k] * Zs
			Im[k] = Vrec[k] / Zm
			Ia[k] = Vrec[k] / Za
			
			theta = math.atan2(Ia[k].imag, Ia[k].real) - math.atan2(Im[k].imag, Im[k].real)
			torque_mec[k] = 2 * a * math.sqrt(Im[k].real**2 + Im[k].imag**2) * math.sqrt(Ia[k].real**2 + Ia[k].imag**2) * (Rb + Rf) * math.sin(theta) / omega_s
			Tload_dy[k] = Tload * (omega_m[k] / (omega_s * (1 - slip_steady)))**2
			omega_m[k + 1] = (dT / J) * (torque_mec[k] - Tload_dy[k]) + omega_m[k]
			
			VIAng[k] = math.atan2(Vrec[k].imag, Vrec[k].real) - math.atan2(I[k].imag, I[k].real)
			PF[k] = math.cos(VIAng[k])
			Pin[k] = math.sqrt(Vrec[k].real**2 + Vrec[k].imag**2) * math.sqrt(I[k].real**2 + I[k].imag**2) * PF[k]
		else:
			# When the rotating speed is large enough, the auxiliary winding is disconnected      
			slip[k] = (omega_s - omega_m[k]) / omega_s
			Zf = 1j * Xem * (R2m / slip[k] + 1j * X2m) / (R2m / slip[k] + 1j * (Xem + X2m))
			Zb = 1j * Xem * (R2m / (2 - slip[k]) + 1j * X2m) / (R2m / (2 - slip[k]) + 1j * (Xem + X2m))
			Rf = Zf.real
			Rb = Zb.real
			
			# Capacitor - start type
			Zm = R1m + 1j * X1m + Zf + Zb
			Za = R1a + 1j * X1a + a**2 * (Zf + Zb) - 1j * Xcr
			Ztotal = Zm * Za / (Zm + Za)
			I[k] = Vs / (Ztotal + Zs)
			Vrec[k] = Vs - I[k] * Zs
			Im[k] = Vrec[k] / Zm
			Ia[k] = Vrec[k] / Za
			
			theta = math.atan2(Ia[k].imag, Ia[k].real) - math.atan2(Im[k].imag, Im[k].real)
			torque_mec[k] = 2 * a * math.sqrt(Im[k].real**2 + Im[k].imag**2) * math.sqrt(Ia[k].real**2 + Ia[k].imag**2) * (Rb + Rf) * math.sin(theta) / omega_s
			Tload_dy[k] = Tload * (omega_m[k] / (omega_s * (1 - slip_steady)))**2
			omega_m[k + 1] = (dT / J) * (torque_mec[k] - Tload_dy[k]) + omega_m[k]
			
			VIAng[k] = math.atan2(Vrec[k].imag, Vrec[k].real) - math.atan2(I[k].imag, I[k].real)
			PF[k] = math.cos(VIAng[k])
			Pin[k] = math.sqrt(Vrec[k].real**2 + Vrec[k].imag**2) * math.sqrt(I[k].real**2 + I[k].imag**2) * PF[k]
	# Transform the dynamic power parameters into equivalent impedance for load flow calculation
	Pstrt = [0] * 3
	PFstrt = [0] * 3
	for j in range(3):
		index = 30 * j + 2 - delay
		if index < 1:
			Pstrt[j] = 1
			PFstrt[j] = PFapp
		else:
			Pstrt[j] = 0.9 * Pin[index - 1]
			PFstrt[j] = PFapp
	return Pstrt, PFstrt

def Zeros(row, col):
	Zeros = []
	for i in range(row):
		col_zeros = []
		for j in range(col):
			col_zeros.append(0.0)
		Zeros.append(col_zeros)
	return Zeros

def OutputFileCSV(filename, Matrix):
	f = open(os.path.dirname(os.path.abspath(__file__)) + r'/results/' + filename,'w')
	for item in Matrix:
		for i in range(len(item)):
			f.write(str(real(item[i])) + ', ' + str(imag(item[i])))
			if i < len(item) - 1:
				f.write(', ')
		f.write('\n')
	f.close()
	return