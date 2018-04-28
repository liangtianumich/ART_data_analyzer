"""
this event selector module filter each event based on various criterias
each criteria is an independent function
"""
import pandas as pd

def df_to_dict(df):
	"""
	this function converts the pandas.Dataframe into a dictionary
	with the key value being the index of pd.Dataframe and the values being
	the values of various columns
	"""
	_final_dict = dict()
	for index,row in df.iterrows():
		_final_dict[index] = row.tolist()
	
	return _final_dict
		
	
def event_select_accept(path_to_file = None):
	"""
	this function filter a single event, i.e. initial configuration, saddle
	configuration and final configuration based on whether this event has been accepted
	or not based on the event_list file inside the ART results
	
	this function should only be used in the exe script file, not to be intertwined
	with the calculator module
	
	Input arguments:
		path_to_file: dir
			this path contains the events.list file in current test
			each test has a different path
	return:
		accepted_events: 
			a dictionary with the key being the id of the attempt
			the values being a list containing the initial, saddle, final
			file string
	for example, it will return 
	{0: ['min1000', 'sad1001', 'min1001', 'accepted'], 
	1: ['min1001', 'sad1002', 'min1002', 'accepted'], 
	...
	
	Note: maybe use an extra input argument event_list_df=None as the 2nd 
	option since event_list_df can be passed as both input and output
	to maintain the consistenency of the code for various functions
	"""
	if path_to_file == None:
		raise Exception("no path of events.list data file has been specified, please specifiy \
		the correct path to the data file")
	events = pd.read_csv(path_to_file,sep='\s+', header=None)
	events.columns = ["ini","sad","fin","status"]
	accepted_events = events.loc[events["status"] == "accepted"]
	
	# drop the column called "status" since now all events are accepted
	# accepted_events = accepted_events.drop('status', axis=1)
	return df_to_dict(accepted_events)

def event_3_criteria():
	df = df.loc[df["status"] == "accepted"]
	df["satisfy_3_criteria"] = False
	
	for index,row in df.iterrows();
		if single_event_3_criteria(row) is True:
			df["satisfy_3_criteria"] = True
	return df		

def single_event_3_criteria(event_ds):
	"""
	maybe organize the following two functions into two methods in a class event
	
	"""
	
	ini, sad, fin = event_list_to_file_path(path_to_data, event_ds)
	
	
def event_list_to_file_path(path_to_data, event_ds):
	"""
	input:
		path_to_data: str
			path string to the data dir
		event_ds: pandas.Series
			single event pandas series containing the str of
			initial configuration in column ini, saddle configuration in column sad,
			final configuration in column fin, 
		
	return:
		ini: pandas.Dataframe
			data of initial configuration
		sad: pandas.Dataframe
			data of saddle configuration
		
		fin: pandas.Dataframe
			data of final configuration
	"""
	
	
	
path_to_file = "~/Documents/alabama/ART_data_analyzer/examples/example_1/events.list"
print event_select_accept(path_to_file)

	
	
