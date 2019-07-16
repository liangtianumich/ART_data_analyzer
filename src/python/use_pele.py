import pele
from pele.landscape import database2graph
from pele.utils.disconnectivity_graph import DisconnectivityGraph
from pele.storage import Database
from data_reader import *
import numpy as np
import matplotlib.pyplot as plt
from util import operation_on_events, data_dir_to_test_dir, event_energy

def make_graph(path_to_data_dir,database):
    # make a graph from the database
    graph = database2graph(database)
    # turn the graph into a disconnectivity graph
    dg = DisconnectivityGraph(graph)
    dg.calculate()
    print "number of minima:", dg.tree_graph.number_of_leaves()
    dg.plot()
    path_to_file = os.path.join(path_to_data_dir,"Disconnectivity_graph.png")
    plt.savefig(path_to_file)
    
def plot_disconnectivity_graph(path_to_data_dir, input_param):
	path_to_db = os.path.join(path_to_data_dir,'pele_database.db')
	db = Database(db=path_to_db, accuracy=0.001, createdb=True)
	
	list_of_test_id = input_param["list_of_test_id"]
	num_of_proc = input_param["num_of_proc"]
	
	operation = lambda x: single_event_adder(x, path_to_data_dir,db)
	
	result_list = operation_on_events(path_to_data_dir, list_of_test_id, operation, num_of_proc = num_of_proc)
	
	make_graph(path_to_data_dir,db)
	
	print ("done adding all final selected events into pele database and plotting their disconnectivity graph!")
	
def single_event_adder(event_state, path_to_data_dir,db):
	if 'test' in event_state[0]:
		test_id = int(event_state[0][4:])
	else:
		test_id = int(event_state[0])
	path_to_test_dir = data_dir_to_test_dir(path_to_data_dir, test_id)
	path_to_curr_result = path_to_test_dir + "/results"
	
	init, sad, fin = event_state[1][0], event_state[1][1], event_state[1][2]
	path_to_curr_event = path_to_curr_result + "/event_" + init + "_" + sad + "_" + fin
	
	path_to_file_ini = path_to_test_dir + '/' + init + ".dump"
	path_to_file_sad = path_to_test_dir + '/' + sad + ".dump"
	path_to_file_fin = path_to_test_dir + '/' + fin + ".dump"

	initial_config_data = read_data_from_file(path_to_file_ini)[['x','y','z']]
	saddle_config_data = read_data_from_file(path_to_file_sad)[['x','y','z']]
	final_config_data = read_data_from_file(path_to_file_fin)[['x','y','z']]
	
	test_eng = event_energy(path_to_test_dir)
	ini_eng = test_eng[init]
	sad_eng = test_eng[sad]
	fin_eng = test_eng[fin]
	
	minimum1 = db.addMinimum(ini_eng, np.array(initial_config_data))
	minimum2 = db.addMinimum(fin_eng, np.array(final_config_data))
	trans = db.addTransitionState(sad_eng, np.array(saddle_config_data), minimum1, minimum2)
	return
