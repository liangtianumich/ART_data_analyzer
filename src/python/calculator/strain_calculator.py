import numpy as np
import pandas as pd
import scipy.spatial as spatial
import os
import pickle
import imp
strain_calculator = imp.load_source('module.name', '../src/python/calculator/periodic_kdtree.py')
PeriodicCKDTree = strain_calculator.PeriodicCKDTree
#from periodic_kdtree import PeriodicCKDTree

class Atom(object):
	"""
	this class initialize the atomic coordinates of single atom for later operation
	the class object can be manipulated as regular numerical operation by using python operator overloading, 
	see http://thepythonguru.com/python-operator-overloading/
	"""
	def __init__(self, atom_loc, box_dim=None, atom_id=1,item=1):
		"""
		constructor intialize the properties of an atom as state variables of the atom, 
		such as atomic coordinates, atom_id, item number in the simulation system
		the size of simulation system
		
		Arguments:
			atom_loc: list or np.array
				list that stores atomic coordinate x,y,z
			box_dim: list, optional, default None object
				list that stores simulation box size lx,ly,lz
		"""
		if type(atom_loc) != list or len(atom_loc) != 3:
			raise Exception("atom location coordinates incorrect")
		if box_dim != None:
			if type(box_dim) != list or len(box_dim) != 3:
				raise Exception("simulation box dimension incorrect")
		self.atom_loc = atom_loc
		self.box_dim = box_dim
		self.atom_id = atom_id
		self.item = item
		self.Z = 
	
	def __add__(self, other):
		final_atom_loc = (np.array(self.atom_loc) + np.array(other.atom_loc)).tolist()		
		return Atom(final_atom_loc)
	
	def __sub__(self, other):
		final_atom_loc = (np.array(self.atom_loc) - np.array(other.atom_loc)).tolist()
		return Atom(final_atom_loc)
	
	def __mul__(self, other):
		final_atom_loc = (np.array(self.atom_loc) * np.array(other.atom_loc)).tolist()
		return Atom(final_atom_loc)
	
	def __rmul__(self, other):
		return self.__mul__(other)
	
	def __div__(self, other):
		final_atom_loc = (np.array(self.atom_loc) / np.array(other.atom_loc)).tolist()
		return Atom(final_atom_loc)
	
	def __rdiv__(self, other):
		return self.__div__(other)
	
	@classmethod
	def distance(cls, atom_1, atom_2):
		"""
		class method to calculate the distance between atom_1 and atom_2,
		where atom_1 and atom_2 are the instances of class Atom
		
		Arguments:
			atom_1: instance of class Atom 
				atomic coordinates of atom_1
			atom_2: instance of class Atom
				atomic coordinates of atom_2
		return:
			distance: float
				the regular distance between atom_1 and atom_2 without considering
				pbc boundary condition
		"""
		
		return np.linalg.norm((atom_1-atom_2).atom_loc)
	
	@classmethod
	def distance_pbc(cls, atom_1, atom_2):
		"""
		class method to calculate the distance between atom_1 and atom_2 under
		periodic boundary condition, where atom_1 and atom_2 are the 
		instances of class Atoms
		
		Arguments:
			atom_1: instance of class Atom 
				atomic coordinates of atom_1
			atom_2: instance of class Atom
				atomic coordinates of atom_2
		return:
			distance: float
				the minimum pair distance between atom_1 and atom_2 under pbc
		"""
		if atom_1.box_dim is None or atom_2.box_dim is None:
			raise Exception("simulation box size has not been specified")
		if atom_1.box_dim != atom_2.box_dim:
			raise Exception("simulation box size does not match")
		
		[lx,ly,lz] = [atom_2.box_dim[0],atom_2.box_dim[1],atom_2.box_dim[2]]
		
		_pair_list = np.array([[0,0,0],[lx,0,0],[-lx,0,0],[0,ly,0],[0,-ly,0],[0,0,lz],[0,0,-lz]])
		
		_pair_distance = []
		
		for _pair in _pair_list:
			_curr_pair_distance = Atom.distance(atom_1, Atom(atom_2.atom_loc + _pair))
			_pair_distance.append(_curr_pair_distance)
		return min(_pair_distance)
	
	@classmethod
	def from_ds(cls, data):
		"""
		create a Atom class object from single column pandas.DataFrame or pandas.Series
		Input argument:
			data: instance of single column pandas.DataFrame or pandas.Series
				data of a single atom
		"""
		if isinstance(data, pd.DataFrame) == False and isinstance(data, pd.Series) == False:
			raise Exception("data converted into Atom class object must be single column pandas.DataFrame or pandas.Series")
		data = data.squeeze()
		_x,_y,_z = data["x"], data["y"], data["z"]
		_atom_id = data["atom_id"]
		_item = data["item"]
		
		return Atom(atom_loc = [_x,_y,_z],atom_id = _atom_id, item = _item)
	
	@classmethod
	def classify_df(cls, data):
		"""
		classify/group all atoms inside configuration data based on their property
		atom_type
		return:
			groups: a dict
				a dictionary with key is atom_id integer, value is the subset of
				pandas.Dataframe of all atoms with this atom_id		
		"""
		if isinstance(data, pd.DataFrame) == False:
			raise Exception("data must be pandas.Dataframe")
		#get unique atom_type id and sorting
		unique_atom_type = sorted(data["atom_id"].unique())
		# find the subset dataframe for each atom_type
		# put their into a dictionary
		# tuple pair key, val in .items() might be useful
		groups = dict()
		for i in unique_atom_type:
			groups[i] = data.loc[data["atom_id"] == i]
		return groups
	
	@classmethod
	def to_list(cls, data):
		"""
		this method convert a class Atom object into a list whose three elements
		are x,y,z
		"""
		if isinstance(data, Atom) == False:
			raise Exception("data must be a class object")
		x,y,z = (data.atom_loc)[0], (data.atom_loc)[1], (data.atom_loc)[2]
		return [x,y,z]


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
	
	

def NN_finder_all(initial_config_data,cut_off_distance, box_dim, atom_list = None, save_results = True):
	"""
	A very general nearest neigbor finder function calculate multiple atom's nearest neighbor all at once using
	the efficient cKDTree algorithm, the multiple atoms whose item number 
	is listed inside the atom_list input argument,
	the default is to calculate all atoms inside initial_config_data file
	User can customize which atoms to calculate by specifying in atom_list
	Input arguments:
	initial_config_data: instance of pandas.Dataframe
		configuration data
	
	cut_off_distance: dict
		dictionary contains the multiple pair cut-off distance
		currently use tuples as keys for immutability, frozenset may be another way
		but it reduce duplicates
		in order to access the turple key without order preference, convert
		
		https://stackoverflow.com/questions/36755714/how-to-ignore-the-order-of-elements-in-a-tuple
		https://www.quora.com/What-advantages-do-tuples-have-over-lists
		For example,
		{(1,1):3.7,(1,2):2.7,(2,2):3.0} means that atom_type 1 and 1 cut-off
		is 3.7, etc
	
	box_dim: list
		a list containing the spatial dimension of simulation box size in x, y, z
		
	atom_list: list
		the list containing the item number of interested atoms whose nearest neighbors
		are being found
	
	save_results: boolean, default True
		specify whether to save the results dictionary into a pkl file
	
	Note:
	this cKDtree algorithm is efficient when:
	you have many points whose neighbors you want to find, you may save 
	substantial amounts of time by putting them in a cKDTree and using query_ball_tree
	
	for molecular simulation: 
	https://github.com/patvarilly/periodic_kdtree
	
	returns:
		nn: dict()
			key is item id of interested atom
			
			values is the pandas.Dataframe of nearest neighbor for atom
			of interest
	"""
	# set up path_to_file and check results out of this function before calling it
	# if check_results is True: 
	# if path_to_file is None or os.path.exists(path_to_file):
	# raise Exception("NN results file not found, please specify the correct path to the file")
		
		
	# if there is no atom_list specified, use all atoms in initial_config_data
	if atom_list is None:
		atom_list = (initial_config_data["item"]).tolist()
	
	_data = initial_config_data
	
	groups = Atom.classify_df(_data)
	
	#_atom_data = initial_config_data[['x','y','z']]
	
	_interested_data = _data.loc[_data['item'].isin(atom_list)]
	
	interested_groups = Atom.classify_df(_interested_data)
	
	#_interested_atom = _interested_data[['x','y','z']]
	
	nn = dict()
	# build the efficient nearest neighbor KDTree algorithm
	# default distance metric Euclidian norm p = 2
	# create tree object using the larger points array
	for (i, int_group) in interested_groups.items():
		for (j, atom_group) in groups.items():
			# comparing atom_type_i and atom_type_j
			for pair in [(i,j),(j,i)]:
				if pair in cut_off_distance:
					 curr_cut_off = cut_off_distance[pair]
			
			# iterate over each row seems inefficient for (index, curr_atom) in int_group.iterrows()
			result_tree = PeriodicCKDTree(box_dim, atom_group[['x','y','z']].values)
			result_groups = result_tree.query_ball_point(int_group[['x','y','z']].values, curr_cut_off)
			#indices = np.unique(IT.chain.from_iterable(result_groups))
			
			#for (int_NN,(index,int_atom)) in (result_groups,int_group.iterrows()):
			k = 0
			for index,int_atom in int_group.iterrows():
				# int_NN is a list of index of NN, index is according to the order
				# in atom_group 
				# curr_NN is a dataframe storing NN found for current atom_group
				int_NN = result_groups[k]
				curr_NN = atom_group.iloc[int_NN]
				if int_atom["item"] not in nn:
					nn[int_atom["item"]] = curr_NN
				elif int_atom["item"] in nn:
					nn[int_atom["item"]] = nn[int_atom["item"]].append(curr_NN)
				k = k + 1
	
	# it is best practice to save this NN dictionary results into a pkl file 
	# to prevent rerun, if this file exists, let user know that
	# the file_of_nearest_neighbor exists before calling it
	if save_results is True:
		with open('nn_result.pkl', 'w') as f:
			pickle.dump(nn,f)
			f.close()
	return nn
	
	
def NN_finder(initial_config_data,atom_item_id,cut_off_distance):
	
	"""
	this function is used to only find the nearest neighbor of a single atom;
	this function take the pandas.Dataframe of initial configuration, returns
	the subset pandas.Dataframe belongs to the nearest neighbors of the atom_id
	
	input arguments:
		initial_data: instance of pandas.Dataframe
			pandas.Dataframe of complete initial atomic configuration
		
		atom_id: integer
			the interested atom item_id that is an integer storing in the pandas.Dataframe item column
			assuming the same item id refers to the same atom after deformation
		
		cut_off_distance: dict
			a dictionary containing the multipe pair cut-off distances in A, 
			might be interested to be extended to triplets in the future etc
			
			Example:
			{"1-2": 2.8; "1-3": 3.2} shows that atom_type_1 and atom_type_2 pair cut-off distance
			is 2.8 A and atom_type 1 and atom_type 3 is 3.2
			
		Note: cut_off_distance can be obtained as a constant/constants approximated 
		for all configurations as the first min distance of the initial quenched configuration 
	
	returns:
		item_id: list of integers
			the item_id of the found nearest neighbours of interested atom_id
		
		nn_data: instance of pandas.Dataframe
			pandas.Dataframe that store all nearest neighbor atoms of 
			interested atom with atom_id
	"""
	df = initial_config_data
	
	# pandas.Series of interested atom
	interested_atom_df = df[df["item"] == atom_item_id]
	
	df["atom_id"].uniques
	interested_atom = Atom.from_ds(interested_atom_df)
	
	# one inefficient brutal force way to code is to iterate over all dataframe and 
	# calculate the distance_pbc and compared based on pair atom_types
	for (_index, _row) in df.iterrows():
		curr_atom = Atom.from_ds(_row)
		#if curr_atom.atom_id ==
	
	
	
	
	
def pdf_calculator_pair(initial_config_data,atom_type_1,atom_type_2):
	"""
	this function calcualte the pair distribution function from scratch in a rigorous way
	it takes the configurational data stored inside the pandas.Dataframe
	initial_config_data, calculate the pair distribution function of pair atom_type_1,
	atom_type_2
	
	Input Arguments:
	initial_config_data: pandas.Dataframe
		configurational data
	atom_type_1: integer
		integer for atom_type_1
	atom_type_2: integer
		integer for atom_type_2
		
	Note: pair distance is calculated as the minimum under periodic boundary condition
	by the distance_pbc method of Atom class
	"""
	atom_1 = pd.DataFrame([])
	atom_2 = pd.DataFrame([])
	
	for (index,atom) in initial_config_data.iterrows():
		if atom["atom_id"] == atom_type_1:
			atom_1 = atom_1.append(atom)
		if atom["atom_id"] == atom_type_2:
			atom_2 = atom_2.append(atom)
	for (in_1, a1) in atom_1.iterrows():
		for (in_2,a2) in atom_2.iterrows():
			at_1 = Atom.from_ds(a1)
			at_2 = Atom.from_ds(a2)
			Atom.distance_pbc(at_1, at_2)
			# to be continued using binning if necessary...
			
		

def first_min_pdf(pdf):
	"""
	this function take the pair distribution function and returns its
	first minumum distance
	"""
	# to be continued in the future if necessary for exact accuracy
	return None
	
def cut_off_distance_calculator(config_data):
	"""
	this function takes a atomic configuration pandas.Dataframe,
	calculate all possible atomic_type pairs' pair distribution function,
	get the cut-off distance for each specific atomic_type pair as the first minumum
	of each pair distribution function,
	
	output the cut-off distances of each pair as a dictionary in the order
	[atom_type_1-atom_type_2,atom_type_1-atom_type_3,atom_type_2-atom_type_3..]
	"""
	atom_type = ((config_data["atom_id"].unique()).astype(np.int64)).tolist()
	atom_type.sort()
	config_size = config_data["item"].size
	cut_off_distance = []
	for i in range(len(atom_type)):
		for j in range(i+1,len(atom_type)-1):
			pair_dist_fun = pdf_calculator_pair(config_data, atom_type[i],atom_type[j])
			cut_off_distance.append(first_min_pdf(pair_dist_fun))
	return cut_off_distance


