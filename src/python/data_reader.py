"""
this data reader module locates and reads the correct ARC data files \
(such as dump files) that stores the atomic coordinates etc
current implementation includes read dump files, more file types can be 
anticipated in the future
"""
import os
import pandas as pd

def read_data_from_file(path_to_file=None, quiet = False):
	try:
		return read_data_from_dump(path_to_file, quiet)
	except IOError:
		path_to_file = path_to_file[:-5]
		return read_data_from_non_dump(path_to_file, quiet)
		
def read_data_from_dump(path_to_file=None, quiet = False):
	"""
	this function takes the full path to the dump data file and read the data file into
	a pandas.Dataframe object for futher data analysis
	
	output: an instance of pandas.Dataframe
		pandas.Dataframe with columes = ['item', 'atom_id', 'x','y','z']
		index = range(#_of_atoms)
	
	"""
	if path_to_file==None:
		raise Exception("no path of data file has been specified, please specifiy \
		the correct path to the data file")
	else:
		data = pd.read_csv(path_to_file,sep='\s+',skiprows = 8)
		result = data.iloc[:,0:5]
		if quiet is False:
			print "read dump data file:", path_to_file
		result.columns = ['item', 'atom_id', 'x','y','z']
		return result

def read_data_from_non_dump(path_to_file=None, quiet = False):
	"""
	this function takes the full path to the non-dump refconfig-type data file and read the data file into
	a pandas.Dataframe object for futher data analysis
	
	output: an instance of pandas.Dataframe
		pandas.Dataframe with columes = ['item', 'atom_id', 'x','y','z']
		index = range(#_of_atoms)
	
	"""
	if path_to_file==None:
		raise Exception("no path of data file has been specified, please specifiy \
		the correct path to the data file")
	else:
		data = pd.read_csv(path_to_file,sep='\s+',skiprows = 2)
		result = data.iloc[:,0:4]
		if quiet is False:
			print "read non-dump data file:", path_to_file
		result.columns = ['atom_id', 'x','y','z']
		result['item'] = result.index +1
		return result

def read_data_from_lammps_data(path_to_file=None, quiet = False):
	"""
	this function takes the full path to the lammps data file output by Ovito and read the data file into
	a pandas.Dataframe object for futher data analysis
	
	output: an instance of pandas.Dataframe
		pandas.Dataframe with columes = ['item', 'atom_id', 'x','y','z']
		index = range(#_of_atoms)
	
	"""
	if path_to_file==None:
		raise Exception("no path of data file has been specified, please specifiy \
		the correct path to the data file")
	else:
		data = pd.read_csv(path_to_file,sep='\s+',skiprows = 9, header=None)
		result = data.iloc[:,0:5]
		if quiet is False:
			print "read lammps data file:", path_to_file
		result.columns = ['item', 'atom_id', 'x','y','z']
		return result

def read_box_range_from_file(path_to_input_files, sample_name, sample_type):
	
	path_to_file = os.path.join(path_to_input_files, sample_name)
	if sample_type == "lammps_data":
		return read_box_range_from_lammps_data(path_to_file)
	elif sample_type == "dump":
		return read_box_range_from_dump(path_to_file)
	else:
		raise Exception("sample type must be either lammps_data or dump")


def read_box_range_from_lammps_data(path_to_file):
	"""
	this function takes the full path to the lammps data file output by Ovito and read the data file into
	a pandas.Dataframe object for futher data analysis
	
	output: a list
		[[x_low,x_high],[y_low, y_high],[z_low,z_high]]
	
	"""
	data = pd.read_csv(path_to_file,sep='\s+', header=None)
	result = (data.iloc[3:6, 0:2]).astype('float')
	return [[result.iloc[0,0],result.iloc[0,1]],[result.iloc[1,0],result.iloc[1,1]],[result.iloc[2,0],result.iloc[2,1]]]

def read_box_range_from_dump(path_to_file):
	"""
	this function takes the full path to the lammps dump file
	
	output: a list
		[[x_low,x_high],[y_low, y_high],[z_low,z_high]]
	"""
	with open(path_to_file) as f:
		lines = f.readlines()
	box_content = [line.strip() for line in lines][5:8]
	l_range = []
	for line in box_content:
		range_list = []
		for char in line.split():
			range_list.append(float(char))
		l_range.append(range_list)
	return l_range
