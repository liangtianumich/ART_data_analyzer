import numpy as np
import pandas as pd
import scipy.spatial as spatial
import os
import json

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
		create a Atomic_Coordinate object from pandas.Series
		Input argument:
			data: instance of pandas.Series
				data of a single atom
		"""
		if isinstance(data, pd.Series) == False:
			raise Exception("data converted into Atom class object must be pandas.Series")
		[_x,_y,_z] = data["x"], data["y"], data["z"]
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
		if isinstance(data, pd.Dataframe) == False:
			raise Exception("data must be pandas.Dataframe")
		#get unique atom_type id and sorting
		unique_atom_type = sorted(data["atom_id"].unique())
		# find the subset dataframe for each atom_type
		# put their into a dictionary
		# tuple pair key, val in .items() might be useful
		groups = dict()
		for i in unique_atom_type:
			groups[i] = data.loc[[data["atom_id"] == i]]
		return groups


def local_strain_calculator(initial_config_data, saddle_config_data):
	

def NN_finder_all(initial_config_data,cut_off_distance, atom_list = None, save_results = True):
	"""
	A very general nearest neigbor finder function calculate multiple atom's nearest neighbor all at once using
	the efficient cKDTree algorithm, the multiple atoms whose item number 
	is listed inside the atom_list input argument，
	the default is to calculate all atoms inside initial_config_data file.
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

	atom_list: list
		the list containing the item number of interested atoms whose nearest neighbors
		are being found
	
	save_results: boolean, default True
		specify whether to save the results dictionary into a json file
	
	Note:
	this cKDtree algorithm is efficient when:
	you have many points whose neighbors you want to find, you may save 
	substantial amounts of time by putting them in a cKDTree and using query_ball_tree
	
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
			result_tree = spatial.cKDTree(atom_group[['x','y','z']].values)
			result_groups = result_tree.query_ball_point(int_group[['x','y','z']].values, curr_cut_off)
			#indices = np.unique(IT.chain.from_iterable(result_groups))
			
			for (int_NN,(index,int_atom)) in (results_groups,int_group.iterrows()):
				# int_NN is a list of index of NN, index is according to the order
				# in atom_group 
				# curr_NN is a dataframe storing NN found for current atom_group
				curr_NN = atom_group.iloc[int_NN]
				if int_atom["item"] not in nn:
					nn[int_atom["item"]] = curr_NN
				elif int_atom["item"] in nn:
					nn[int_atom["item"]] = nn[int_atom["item"]].append(curr_NN)
	
	# it is best practice to save this NN dictionary results into a json file 
	# to prevent rerun, if this file exists, let user know that
	# the file_of_nearest_neighbor exists before calling it
	if save_results is True:
		with open('nn_result.json', 'w') as f:
			json.dump(nn,f)
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
		if curr_atom.atom_id ==
	
	
	
	
	
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

def local_strain_calculator_atom(NN_data_initial, NN_data_saddle, atom_id):
	"""
	this function calculate the local atomic strain according to the algorithm
	proposed in http://li.mit.edu/Archive//Graphics/A/annotate_atomic_strain/Doc/main.pdf
	
	Input arguments:
		
		initial_data: an instance of pandas.Dataframe
			intial configuration stores atom id, atom type and x,y,z coordinates of the dump file
		saddle_data: an instance of pandas.Dataframe
			saddle configuration stores atom id, atom type and x,y,z coordinates of the dump file

	return:
		a list of pandas.Series instances
		stored various local atomic strain in pandas.Series
		such as the von mises strain, 
	
	explain the details of this funtion, what is the input and output, is the output \
	strain matrix relative to a single atom or for all atoms  */
	"""
