import numpy as np
import pandas as pd
import os
import pickle
from util import Atom, NN_finder_all

def local_strain_calculator_orth(initial_config_data, saddle_config_data, cut_off_distance, box_dim, atom_list = None, save_results = True):
	"""
	this function calculate various local atomic strain quantities for atoms whose item id stored inside
	the atom_list under periodic boundary condition for orthogonal simulation box
	input arguments:
	
	initial_config_data: instance of pandas.Dataframe
		intial configuration data
	
	saddle_config_data: instance of pandas.Dataframe
		saddle configuration data
		
	cut_off_distance: dict
		dictionary contains the multiple pair cut-off distance
		currently use tuples as keys for immutability
		for example, {(1,1):3.7,(1,2):2.7,(2,2):3.0} means that atom_type 1 and 1 cut-off
		is 3.7, etc
		It is suggested that cut_off_distance are used as a user input,
		not implemented by calculating from the first min of pair distribution function
		of initial_config_data since it is not a module that reused a lot
		it is good enough to get this single value through Ovide coordination analysis rdf or other software
		instead of a rigourous calculation
	
	box_dim: list,
		the spatial dimension of simulation box in [x,y,z]
	
	atom_list: list,
		a list storing the atom item id of interested atoms
		default calculating all atoms in the initial_config_data
	
	save_results: boolean, default True
		if True, save calculated local strains results into a pkl file
		
	returns:
		strain: dict()
			a dictionary with key being the item id number of the interested atom
			values being the calculated various local atomic strains measures 
			stored inside in a list
			
	Note:
		the nearest neighbor is determined through the atomic configuration in
		initial_config_data
	
	
	"""
	# if not specifying the atom_list, choose all atoms in the initial_config_data
	if atom_list is None:
		atom_list = (initial_config_data["item"]).tolist()
	_data = initial_config_data
	
	# check if the nn results.pkl file exists or not
	
	nn = NN_finder_all(initial_config_data, cut_off_distance, box_dim, atom_list)
	
	strain = dict()
	
	for item in atom_list:
		
		item_nn = nn[item]
		NN_list = item_nn["item"].tolist()
		NN_list.append(item)
		
		# NN_intial is the pandas.Dataframe of intial data with interested_atom and its NN
		NN_initial = initial_config_data.loc[initial_config_data["item"].isin(NN_list)]
		
		# NN_saddle is the pandas.Dataframe of saddle data with interested_atom and its NN
		NN_saddle = saddle_config_data.loc[saddle_config_data["item"].isin(NN_list)]
		#NN_intial = initial_config_data.loc[(initial_config_data['item'].isin(item_nn["item"].tolist()))]
		#NN_intial.append(item)
		#NN_saddle = (saddle_config_data['item'].isin(item_nn["item"])).tolist() 
		#saddle_config_nn = saddle_config_data.loc[(saddle_config_data['item'].isin(item_nn["item"])).tolist()]
		
		# local_strains should be a dict as well since it has multiple output strain
		# or a list
		local_strains = local_strain_calculator_atom_orth(NN_initial, NN_saddle, item, box_dim)
		strain[item] = local_strains
	
	if save_results is True:
		with open('strain_result.pkl', 'w') as f:
			pickle.dump(strain,f)
			f.close()
	return strain


def local_strain_calculator_atom_orth(initial_config_atom, saddle_config_atom, atom_item, box_dim):
	"""
	this function calculate the local atomic strain according to the algorithm
	proposed in http://li.mit.edu/Archive//Graphics/A/annotate_atomic_strain/Doc/main.pdf
	
	Input arguments:
		
		initial_config_atom: an instance of pandas.Dataframe
			interested atom and its NN atoms in intial configuration data stores atom id, atom type and x,y,z coordinates of the dump file
		
		saddle_config_atom: an instance of pandas.Dataframe
			interested atom and its NN atoms in saddle configuration data stores atom id, atom type and x,y,z coordinates of the dump file
		
		atom_item: integer
			the item id of interested atom
			
		box_dim: list
			the simulation box spatial dimension to implement periodic boundary condition
	return:
		atom_strain: list
			storing various local atomic strains in list
			such as the von mises shear strain invariants, hydrostatic invariant 
	"""
	# calculate dij and d0_ij
	NN_initial = initial_config_atom.loc[initial_config_atom['item'] != atom_item]
	Atom_initial = initial_config_atom.loc[initial_config_atom['item'] == atom_item]
	Atom_ini_obj = Atom.from_ds(Atom_initial)
	
	
	NN_saddle = saddle_config_atom.loc[saddle_config_atom['item'] != atom_item]
	Atom_saddle = saddle_config_atom.loc[saddle_config_atom['item'] == atom_item]
	Atom_sad_obj = Atom.from_ds(Atom_saddle)
	
	Dim = 3
	V = np.zeros(shape=(Dim,Dim))
	W = np.zeros(shape=(Dim,Dim))
	
	for (index,atom_ini_NN) in NN_initial.iterrows():
		# d0_ji in pandas.Series
		d0_ji = Atom.from_ds(atom_ini_NN) - Atom_ini_obj
		d0_ji = Atom.to_list(d0_ji)
		atom_sad_NN = Atom.from_ds(NN_saddle.loc[NN_saddle["item"] == atom_ini_NN["item"]])
		# d_ji in pandas.Series
		d_ji = atom_sad_NN - Atom_sad_obj
		d_ji = Atom.to_list(d_ji)
		
		# begin implement pbc for d_ji and d0_ji as in 
		# https://en.wikipedia.org/wiki/Periodic_boundary_conditions
		# dx = x[j] - x[i];
		# dx -= x_size * nearbyint(dx * x_rsize)
		
		Near_int_d = [int(round(i)) for i in np.array(d_ji) * 1.0/np.array(box_dim)]
		d_ji = np.array(d_ji) - np.array(box_dim) * np.array(Near_int_d)
		
		Near_int_d0 = [int(round(i)) for i in np.array(d0_ji) * 1.0/np.array(box_dim)]
		d0_ji = np.array(d0_ji) - np.array(box_dim) * Near_int_d0
		
		
		# begin calculate the V and M matrix
		V[0][0] = V[0][0] + d0_ji[0] * d0_ji[0]
		V[0][1] = V[0][1] + d0_ji[0] * d0_ji[1]
		V[0][2] = V[0][2] + d0_ji[0] * d0_ji[2]
		V[1][0] = V[1][0] + d0_ji[1] * d0_ji[0]
		V[1][1] = V[1][1] + d0_ji[1] * d0_ji[1]
		V[1][2] = V[1][2] + d0_ji[1] * d0_ji[2]
		V[2][0] = V[2][0] + d0_ji[2] * d0_ji[0]
		V[2][1] = V[2][1] + d0_ji[2] * d0_ji[1]
		V[2][2] = V[2][2] + d0_ji[2] * d0_ji[2]

		W[0][0] = W[0][0] + d0_ji[0] * d_ji[0]
		W[0][1] = W[0][1] + d0_ji[0] * d_ji[1]
		W[0][2] = W[0][2] + d0_ji[0] * d_ji[2]
		W[1][0] = W[1][0] + d0_ji[1] * d_ji[0]
		W[1][1] = W[1][1] + d0_ji[1] * d_ji[1]
		W[1][2] = W[1][2] + d0_ji[1] * d_ji[2]
		W[2][0] = W[2][0] + d0_ji[2] * d_ji[0]
		W[2][1] = W[2][1] + d0_ji[2] * d_ji[1]
		W[2][2] = W[2][2] + d0_ji[2] * d_ji[2]
	
	J = np.dot(np.linalg.inv(V), W)
	
	mu = (np.dot(J, np.transpose(J)) - np.identity(3)) * 0.5
	
	mu_hydro = np.trace(mu)/3.0
	
	mu_off_diag = (mu - mu_hydro * np.identity(3))
	
	mu_Mises = (0.5 * np.trace(np.dot(mu_off_diag, mu_off_diag))) ** 0.5
	
	return [mu_hydro, mu_Mises]
