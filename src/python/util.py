"""
utilities module to store useful classes and functions for other modules
"""
import numpy as np
import pandas as pd
import pickle
import os
import json
import re
import operator
from pathos.multiprocessing import ProcessingPool as Pool
from data_reader import *
from periodic_kdtree import PeriodicCKDTree
from functools import wraps
import time

def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print ("Total time running %s: %s seconds" %
               (function.func_name, str(t1-t0))
               )
        return result
    return function_timer

def data_dir_to_test_dir(path_to_data_dir, test_id):
	"""
	this function return the path to test dir for test with test_id
	"""
	path_to_curr_test = [path_to_data_dir + "test%s"%test_id, path_to_data_dir + "%s"%test_id]
	if os.path.exists(path_to_curr_test[0]):
		return path_to_curr_test[0]
	elif os.path.exists(path_to_curr_test[1]):
		return path_to_curr_test[1]
	else:
		raise Exception("each test dir in %s should be named either test_1 or 1 for test_id = 1"%path_to_data_dir)
	
	
class Atom(object):
	"""
	this class initialize the atomic coordinates of single atom for later operation
	the class object can be manipulated as regular numerical operation by using python operator overloading, 
	see http://thepythonguru.com/python-operator-overloading/
	"""
	def __init__(self, atom_loc, box_dim=None, atom_id=1,item=1, atom_z=None):
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
		self.atom_z = atom_z
	
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
			_curr_pair_distance = Atom.distance(atom_1, Atom((np.array(atom_2.atom_loc) + _pair).tolist()))
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


class Configuration(object):
	def __init__(self,path_to_config_dump_data,box_dim, quiet = True):
		"""
		
		"""
		
		self.data = read_data_from_file(path_to_config_dump_data, quiet)
		self.box_dim = box_dim
	
	@classmethod
	def distance(cls,config_1, config_2):
		"""
		this method implement the distance in the following definition
		sqrt(sum(dx**2+dy**2+dz**2) for each atom in config)
		Input arguments:
			config_1: instance of Configuration object
				its data attribute returns the pandas.DataFrame
				
			config_2: instance of Configuration object
				its data attribute returns the pandas.DataFrame	
		return: 
			distance: scalar float
				the distance between two atomic configurations
			
		"""
		sorted_data_1 = (config_1.data).sort_values('item')
		sorted_data_2 = (config_2.data).sort_values('item')
		dr = sorted_data_1 - sorted_data_2
		return np.linalg.norm(dr)
	
	@classmethod
	def distance_pbc(cls, config_1, config_2):
		"""
		this method implement the distance in the following definition
		sqrt(sum(dx**2+dy**2+dz**2) for each atom in config)
		the dx, dy, dz should follow the minumum image convention in pbc
		
		Input arguments:
			config_1: instance of Configuration object
				its data attribute returns the pandas.DataFrame
				
			config_2: instance of Configuration object
				its data attribute returns the pandas.DataFrame	
		return: 
			distance: scalar float
				the distance between two atomic configurations under pbc
		
		Note: each atom follow the minimum image convention
		"""
		sorted_data_1 = (config_1.data).sort_values('item')
		sorted_data_2 = (config_2.data).sort_values('item')
		box_dim_1 = config_1.box_dim
		box_dim_2 = config_2.box_dim
		i = 0
		total_distance = 0.0
		for index, row in sorted_data_1.iterrows():
			atom_1 = Atom.from_ds(row)
			atom_1.box_dim = box_dim_1
			atom_2 = Atom.from_ds(sorted_data_2.iloc[i])
			atom_2.box_dim = box_dim_2
			atom_dist = Atom.distance_pbc(atom_1,atom_2)
			total_distance = total_distance + atom_dist ** 2
			i = i + 1
		return total_distance ** 0.5	

def event_distance(path_to_test_dir, event_state, box_dim):
	"""
	this function calculates the distance between final and initial configuration
	of a single event
	"""
	
	event_init,event_sad,event_fin = event_state[0],event_state[1],event_state[2]
	
	init_path = path_to_test_dir + '/' + event_init +".dump"
	sad_path = path_to_test_dir + '/' + event_sad +".dump"
	fin_path = path_to_test_dir + '/' + event_fin +".dump"
	
	init_config = Configuration(init_path,box_dim)
	fin_config = Configuration(fin_path,box_dim)
	
	return Configuration.distance_pbc(init_config,fin_config)

def get_list_of_selected_events_str(path_to_data_dir, list_of_test_id):
	"""
	this function returns a list containing the strings of the selected events
	after stage I two criteria for all tests in list_of_test_id
	"""
	all_selected_events = []
	for i in list_of_test_id:
		path_to_curr_test = data_dir_to_test_dir(path_to_data_dir,i)
		path_to_selected_events = path_to_curr_test + "/results/selected_events.json"
		# skip the test who do not have selected_events.json in all tests specified in
		# list_of_test_id
		if os.path.exists(path_to_selected_events):
			selected_events = json.load(open(path_to_selected_events,'r'))
			# value = [init_state, sad_state, fin_state]
			selected_list = selected_events.values()
			for event in selected_list:
				event_str = ("test%s"%i, [event[0],event[1],event[2]])
				all_selected_events.append(event_str)
	return all_selected_events

def get_list_of_final_filtered_events_str(path_to_data_dir):
	"""
	this function returns a list containing the strings of the selected events
	after both stage I and stage II three criteria for all tests in list_of_test_id
	"""
	all_selected_events = []
	path_to_final_selected_events = path_to_data_dir + "final_selected_events.json"
	if os.path.exists(path_to_final_selected_events):
		final_selected_events = json.load(open(path_to_final_selected_events,'r'))
		# value = [init_state, sad_state, fin_state]
		for event in final_selected_events:
			event_str = (event[0], event[1])
			all_selected_events.append(event_str)
	else:
		raise Exception("events in current list_of_test_id has not been filtered yet! run art_data -s SETTINGS --filter first")
	return all_selected_events
	
def event_act_relax_energy(event, path_to_data_dir):
	"""
	this function calculates the activiation energy for a single event
	from the test/log.file.1
	
	Input:
		path_to_data_dir:
		
		event: a list
		a list with the 1st element being the test_id, e.g test1
		with the 2nd element being a list of event_str,[init, sad, fin]
	
	return:
		a tuple of event_state, and its activation energy
	"""
	if 'test' in event[0]:
		test_id = int(event[0][4:])
	else:
		test_id = int(event[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	
	#path_to_test_dir = path_to_data_dir + event[0]
	test_energy_list = event_energy(path_to_test_dir)
	init_eng = test_energy_list[event[1][0]]
	sad_eng = test_energy_list[event[1][1]]
	fin_eng = test_energy_list[event[1][2]]
	act_eng = sad_eng - init_eng
	relax_eng = sad_eng - fin_eng
	return (event, act_eng, relax_eng)

def state_energy_barrier(path_to_test_dir, init, state):
	"""
	this function calculates the activiation energy for a initial state and 
	another interested state in a test
	from the test/log.file.1
	
	Input:
		path_to_test_dir
		
		init: a str
		
		state: a str
	
	return:
		a tuple of event_state, and its activation energy
	"""
	test_energy_list = event_energy(path_to_test_dir)
	init_eng = test_energy_list[init]
	state_eng = test_energy_list[state]
	eng_barrier = state_eng - init_eng
	return eng_barrier

def event_energy(path_to_test_dir):
	"""
	Input:
		path_to_test_dir: str
			path of log.file.1 for test1 that contains the energy information
	
	Output:
		energy_of_events: a dict
		a dict giving the energy of each state of each successful event
		key being the state str, such as min1000,sad1001 etc
		value being the float value of the energy
	"""
	f=open(path_to_test_dir+"/log.file.1")
	string=f.read()
	energy_of_events = dict()
	pattern = "Configuration[\s]+stored[\s]+in[\s]+file[\s]+:[\s]+([minsad\d]+)\s+.+\s+.+\s+.+\s+.+Starting\s+from\s+minconf.+\n.+Reference\s+Energy\s+.eV.+:\s+([-+.eE\d]+)"
	match = re.search(pattern, string)
	energy_of_events[match.group(1)] = float(match.group(2))
	
	pattern_other = "Configuration[\s]+stored[\s]+in[\s]+file[\s]+:[\s]+([minsad\d]+)\s+.+Total\s+energy\s+[a-zA-Z]+\s+.eV.+:\s+([-+.eE\d]+)"
	match_other = re.findall(pattern_other, string)
	for state,energy in match_other:
		energy_of_events[state] = float(energy)
	return energy_of_events

def event_strain_disp(event_strain_dict,event_disp_dict):
	"""
	this function takes a strain dictionary with the key being atom item_id
	and value being a list of [volume, von_Mises] converts it into a list
	"""
	vol_strain = []
	shear_strain = []
	disp = []
	
	keys = list(event_strain_dict.keys())
	keys.sort()
	for i in keys:
		vol_strain.append(event_strain_dict[i][0])
		shear_strain.append(event_strain_dict[i][1])
		disp.append(event_disp_dict[i])
	return (vol_strain, shear_strain, disp)

def run_tests_triggered_atom_is_max_disp(path_to_data_dir, input_param):
	list_of_test_id = input_param["list_of_test_id"]
	num_of_proc = input_param["num_of_proc"]
	operation = lambda x: triggered_atom_is_max_disp(path_to_data_dir, x)
	results_list = operation_on_events(path_to_data_dir, list_of_test_id, operation, num_of_proc=num_of_proc)
	print "done checking triggered atom is max displacement atom for interested tests!"
	
def triggered_atom_is_max_disp(path_to_data_dir, event):
	"""
	this function check if a single event the maximum displacement aton is
	the triggered atom stored inside the bart.sh
	"""
	print "event:", event[0],event[1][0],event[1][1],event[1][2]
	
	if 'test' in event[0]:
		test_id = int(event[0][4:])
	else:
		test_id = int(event[0])
	
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	
	#path_to_test_dir = path_to_data_dir + event[0]
	
	triggered_atom_index = read_from_art_input_file(path_to_test_dir)
	
	init, sad, fin = event[1][0],event[1][1],event[1][2]
	
	path_to_event = path_to_test_dir + "/results/event_" + init + "_" + sad + "_" + fin
	
	path_to_init_sad = path_to_event + "/init_sad"
	
	path_to_displacement = path_to_init_sad + "/displacement_results_dict.pkl"
	
	if os.path.exists(path_to_displacement):
		print "path to displacement:", path_to_displacement
		event_disp = pickle.load(open(path_to_displacement,'r'))
		index_max_disp = max(event_disp.iteritems(), key=operator.itemgetter(1))[0]
		print "max displacement atom index:", index_max_disp
		if len(triggered_atom_index) == 1:
			if triggered_atom_index[0] == index_max_disp:
				print "True"
				return True
			else:
				print "False"
				return False
		else:
			print "multiple triggering atoms exists!"	
	else:
		print("no displacement data has been calculated in current event")
	
	
def operation_on_events(path_to_data_dir, list_of_test_id, operation, num_of_proc=1):
	"""
	this function perform an operation function on each events listed in tests in
	list_of_test_id
	Input arguments:
		operation: function
			an operation acting on a single event
	"""
	test_id = ["test%s"%i for i in list_of_test_id]
	pool = Pool(processes = num_of_proc)
	path_to_final_selected_events = path_to_data_dir + "final_selected_events.json"
	if os.path.exists(path_to_final_selected_events):
		final_selected_events = json.load(open(path_to_final_selected_events,"r"))
		final_interested_events = []
		for event in final_selected_events:
			if event[0] in test_id:
				final_interested_events.append(event)
	else:
		final_interested_events = []
		for test in list_of_test_id:
			path_to_curr_test = data_dir_to_test_dir(path_to_data_dir, test)
			path_to_test_result = path_to_curr_test +"/results"
			path_to_event_list = path_to_test_result + "/selected_events.json"
			if os.path.exists(path_to_event_list):
				event_list = json.load(open(path_to_event_list,"r"))
				for value in event_list.values():
					event = ["test%s"%test,[value[0],value[1],value[2]]]
					final_interested_events.append(event)
			else:
				print "skip current test:", "test%s"%test, "there is no selected events"
	
	# if function operation has no return value, it will return a list of Nones
	result_list = pool.map(operation,final_interested_events)
	return result_list


def event_local_atom_index(initial_config_data, triggered_atom_list, num_of_involved_atom, path_to_init_sad, box_dim, save_results = True, re_calc = False):
	"""
	this function get the local atom atom_ids as a list of lists (that will be flattern into a single list)
	from the triggered atoms (whose item_ids in triggered_atom_list) item_ids
	and num_of_involved atoms NN
	For example, each trigger atom has k th NN. Then for x triggered atoms,
	it would be a x element list of lists (each list with k elements). 
	It will be flattern into a list with k*N, that will be returned 
	
	As in NN_finder_all, PeriodicCkDtree can be easily extended to
	find kth NN for each of the triggered atoms.
	
	Currently, only one atom is triggered,
	In this case, triggered_atom_list is a single element list
	
	return: 
	triggered_atoms_NN_list: a list 
		comes from a flatten numpy array shape (len(triggered_atom_list)*k, )
		a numpy contains all NN of all triggered atoms preserving the order
		of the appearance of triggered atom in triggered_atom_list
		e.g. 1st k elements in triggered_atoms_NN_list belongs to the 1st
		element of triggered_atom_list, etc
	"""
	
	path_to_local_nn_results = path_to_init_sad + "/local_nn_results_dict.pkl"
	print "the triggering atoms are:", triggered_atom_list
	if re_calc is False:
		if os.path.exists(path_to_local_nn_results):
			print "local nn results already calculated and saved in local_nn_results_dict.pkl file, skip calculation"
			local_nn = pickle.load(open(path_to_local_nn_results,'r'))
			return (np.array(local_nn.values()).flatten()).tolist()
	local_nn = dict()
	
	if triggered_atom_list is None:
		raise Exception("try to calculate the NN for triggered atoms, but no triggered atom has been specified")
	
	_data = initial_config_data
	
	_interested_data = _data.loc[_data['item'].isin(triggered_atom_list)]
	
	result_tree = PeriodicCKDTree(box_dim, _data[['x','y','z']].values)
	# the query method calculated NN includes the triggered atom itself as the 1st NN with distance 0
	distances,locations = result_tree.query(_interested_data[['x','y','z']].values, num_of_involved_atom)
	
	if len(triggered_atom_list) > 1 and num_of_involved_atom >1:
		# the 1st element in the return numpy array with shape (k*len(triggered_atom_list),) 
		# is the triggered atom since the distance to itself is 0, which is the minumum
		# locations are ordered in terms of their increasing distance to the triggered atom
		k=0
		for index,row in _interested_data.iterrows():
			NN_array = np.array((_data.iloc[locations[k]])['item'])
			local_nn[row['item']]= (NN_array).tolist()
			k=k+1
		loc_index = locations.flatten()
		final_locations = np.array((_data.iloc[loc_index])['item']).tolist()
	elif len(triggered_atom_list) == 1:
		
		if type(locations) == int or type(locations) == float:
			loc_index = np.array([locations])
		else:
			loc_index = locations.flatten()
		locations = np.array((_data.iloc[loc_index])['item']).tolist()
		local_nn[triggered_atom_list[0]] = locations
		final_locations = locations
	else:
		# len(triggered_atom_list) >1 and num_of_involved_atom == 1:
		for x in triggered_atom_list:
			local_nn[x] = [x]
		final_locations = triggered_atom_list
	if save_results is True:
		with open(path_to_local_nn_results, 'w') as f:
			pickle.dump(local_nn,f)
			f.close()
	return final_locations

def read_from_art_input_file(path_to_test_dir):
	"""
	this function get the triggered atom index from the art input file bart.sh
	"""
	path_to_test_bart = path_to_test_dir + "/bart.sh"
	f=open(path_to_test_bart)
	string=f.read()
	pattern = "([#\n])setenv[\s]+Central_Atom[\s]+([1-9]+)"
	#pattern  #setenv Central_Atom        1
	match = re.search(pattern, string)
	if match.group(1) == "#":
		return [1]
	elif match.group(1) == "\n":
		return [int(match.group(2))]
	
	
	
	
class results(object):
	def __init__(self, path_test_dir, path_to_data):
		
		# events = event_selector(path_to_data_dir)
		# save events into selected_event_list.json
		
		# for event in events:
		# 	path_test_dir = path_to_data_dir + "event_%s_%s"%event[0],event[1]
		
		self.path_to_test_dir = path_to_test_dir 
		self.path_to_input = self.path_to_test_dir + "/input.json"
		self.path_to_results = self.path_to_test_dir + "/results.pkl"
		
		
		self.data = pickle.load(open(self.path_to_results,'r'))
	
	def statistics(self):
		
		# get the statistics of the various quantities calculated in
		# save the file into a file called results_statistics.json
		
		self.path_to_file_ini = self.path_to_test_dir + "min1000.dump"
		self.path_to_file_sad = self.path_to_test_dir + "min1001.dump"

		self.initial_config_data = read_data_from_file(self.path_to_file_ini)
		self.saddle_config_data = read_data_from_file(self.path_to_file_sad)
		
		Event(self.initial_config_data, self.saddle_config_data)



def NN_finder_all(initial_config_data, cut_off_distance, box_dim, path_to_test_dir, atom_list = None, save_results = False, re_calc = False):
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
	
	path_to_test_dir: str
		str of current test result dir, under it, it save data into nn_results.pkl
		
	atom_list: list
		the list containing the item number of interested atoms whose nearest neighbors
		are being found
	
	save_results: boolean, default True
		specify whether to save the results dictionary into a nn_results_dict.pkl file
	
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
		
	path_to_nn_results = path_to_test_dir + "/nn_results_dict.pkl"
	
	if re_calc is False:
		if os.path.exists(path_to_nn_results):
			print "nn results dictionary already calculated and saved in pkl file, skip calculation"
			return pickle.load(open(path_to_nn_results,'r'))
	nn = dict()
		
	# if there is no atom_list specified, use all atoms in initial_config_data
	if atom_list is None:
		atom_list = (initial_config_data["item"]).tolist()
	
	_data = initial_config_data
	
	groups = Atom.classify_df(_data)
	
	#_atom_data = initial_config_data[['x','y','z']]
	
	_interested_data = _data.loc[_data['item'].isin(atom_list)]
	
	interested_groups = Atom.classify_df(_interested_data)
	
	#_interested_atom = _interested_data[['x','y','z']]
	
	
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
		with open(path_to_nn_results, 'w') as f:
			pickle.dump(nn,f)
			f.close()
	return nn


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
