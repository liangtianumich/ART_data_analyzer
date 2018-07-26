"""
this data reader module locates and reads the correct ARC data files \
(such as dump files) that stores the atomic coordinates etc
current implementation includes read dump files, more file types can be 
anticipated in the future
"""
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
		data = pd.read_csv(path_to_file,sep='\s+',skiprows = 2)
		result = data.iloc[:,0:4]
		if quiet is False:
			print "read non-dump data file:", path_to_file
		result.columns = ['atom_id', 'x','y','z']
		result['item'] = result.index +1
		return result
