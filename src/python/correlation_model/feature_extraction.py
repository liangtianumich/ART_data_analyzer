"""
this feature extraction module design well-crafted, physically meaningful features from atomic coordiantes of
different atomistic structures for the purpose of dimension reduction to avoid curse of dimensionality
"""
import numpy as np
import pickle
from scipy.special import sph_harm
from util import Atom
from util import NN_finder_all

def SOAP_matrix_calculator(configuration_data, cut_off_distance, box_dim, path_to_results = None):
	"""
	this Smooth Overlap of Atomic Positions (SOAP) descriptor matrix calculator 
	takes the atomic coordinates data and atom species data and output a SOAP matrix 
	see reference in Conrad W. Rosenbrock, Eric R. Homer, Gábor Csányi 
	and Gus L. W. Hart, "Discovering the building blocks of atomic systems using
	machine learning: application to grain boundaries"
	
	In results.pkl file, under each atom_id key, add a SOAP_vector key
	it is easy to assemble SOAP_matrix from SOAP vector, no need to
	add a global key. So only use this function when need to calculate
	the SOAP matrix from scratch. Pre-check if ALL atoms's SOAP_vector key
	exists by checking shape.
	
	The atomic density function is a g
	
	"""
	if path_to_results is None:
		path_to_results = "results.pkl"
	
	SOAP_matrix = []
	
	# if the SOAP vector is calculated for all atoms, then assemble the SOAP vector
	# from results.pkl file, move this outside of this function to prevent intertwinned
	if os.path.exists(path_to_results):
		results = pickle.load(open(path_to_results,"r"))
		for atom_id,atom_result in results.items():
			if "SOAP_vector" in atom_result:
				curr_SOAP_vector = atom_result["SOAP_vector"]
				SOAP_matrix.append(curr_SOAP_vector)	
			else:
				raise Exception("SOAP vector has not been calculated for atom item_id %d"%atom_id)
		SOAP_matrix = np.matrix(SOAP_matrix)
		
	
	
	# calculate SOAP matrix from scratch
	for index,row in configuration_data.iterrows():
		curr_atom = Atom.from_ds(row)
		if not os.path.exists(path_to_results):
			curr_NN = NN_finder_all(configuration_data,cut_off_distance, box_dim, curr_atom.atom_id)
		else:
			try:
				curr_NN = results[curr_atom.atom_id]["NN"]
			except KeyError:
				print"atom in results.pkl file has no NN results for atom_id %d"%curr_atom.atom_id
				print "now calculating atom item id %d NN"%curr_atom.atom_id
				# this function default save the NN results into key "NN"
				curr_NN = NN_finder_all(configuration_data,cut_off_distance, box_dim, curr_atom.atom_id)
				
				
		#calculate SOAP vector
		curr_atom_SOAP_vector = SOAP_vector(curr_atom, curr_NN)
		SOAP_matrix.append(curr_atom_SOAP_vector)
		
	return np.matrix(SOAP_matrix)
		

def SOAP_vector(curr_atom, curr_NN, save_results = True):
	"""
	This function calculate the SOAP vector for a single atom whose NN
	are specified, this SOAP vector performs like the property of an atom
	so that it can be saved into the results.pkl file along with the strain
	calculation results and NN results, merge strain results and NN results
	'nn_result.pkl' into one file called results.pkl
	
	Input: 
	curr_atom: pandas.Series
	
	curr_NN: pandas.Dataframe
	
	
	output: list
	  SOAP vector for curr_atom
	
	Note:
	Conrad W. Rosenbrock, Eric R. Homer, Gábor Csányi 
	and Gus L. W. Hart, "Discovering the building blocks of atomic systems using
	machine learning: application to grain boundaries"
	
	1) place a Gaussian particle density function at the location of each atom
	within a local environment sphere around the atom, cut_off 5 A
	2) the total density function produced by the neighbors is projected 
	into a spectral basis consisting of radial basis functions and the spherical harmonics,
	each basis function produces a single coefficient pi in the SOAP vector p for the atom	
	the coefficient is cnlm, coming from the summation of C_nlm from each neighboring atom
	The basis function can be seen in eq 37
	Albert P. Bartók, Risi Kondor, Gábor Csányi
	On representing chemical environments
	"""
	
	
	# assigning each atom_type with a weight factor to scale up the atomic number
	# effect
	curr_atom = Atom.from_ds(curr_atom)
	#calculate the dr
	for index, NN_atom in curr_NN.iterrows():
		# for radial and spherical harmonics basis function
		dr = (NN_atom - curr_atom).atom_loc
		# determine coefficient for basis function
		
		
		# add atom_z attribute into the Atom class
 		atom_type_1 = curr_atom.atom_id
		atom_type_2 = NN_atom.atom_id
		if atom_type_1 == 1:
			# type 1 is Cu, atomic number 29; type 2 is Zr, atomic number 40
		elif atom_type_2 == 2:
			
			
		gaussian_spectral_decomposition(dr, curr_atom.atom_z, NN_atom.atom_z)
		
		
		
		
	
	
