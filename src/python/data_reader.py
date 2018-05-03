"""
this data reader module locates and reads the correct ARC data files \
(such as dump files) that stores the atomic coordinates etc
current implementation includes read dump files, more file types can be 
anticipated in the future
"""
import pandas as pd

def read_data_from_dump(path_to_file=None):
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
		#print "data.head", data.head()
		#print "data.iloc[:,0:5]", data.iloc[:,0:5]
		result = data.iloc[:,0:5]
		result.columns = ['item', 'atom_id', 'x','y','z']
		#.rename(columns=['item', 'atom id', 'x','y','z'])
		#result = pd.DataFrame(data.iloc[:,0:5], columns=['item', 'atom id', 'x','y','z'])
		return result
