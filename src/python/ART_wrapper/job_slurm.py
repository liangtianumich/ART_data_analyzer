import os
import subprocess
import time
import multiprocessing as mp

def run_bart_node(path_to_data_dir, list_of_central_atoms):
	"""
	this function will be used to submit a job to a compute node. This job will
	invoke the parallized art running using all cores of this compute node
	"""
	
	list_of_run_dir = []
	for central_atom in list_of_central_atoms:
		path_to_run_dir = os.path.join(path_to_data_dir, str(central_atom))
		list_of_run_dir.append(path_to_run_dir)
	pool = mp.Pool(processes=mp.cpu_count())
	pool.map(run_bart_core,list_of_run_dir)

def run_bart_core(path_to_run_dir):
	os.chdir(path_to_run_dir)
	subprocess.check_call("chmod +x mod_bart.sh; ./mod_bart.sh", shell=True)
	print "done runing art for %s"%path_to_run_dir

central_atoms = os.environ["ATOM_LIST"][1:-1].split(" ")
list_of_central_atoms = map(int, central_atoms)
path_to_data_dir = os.environ["DATA_DIR"]

start_time = time.time()
run_bart_node(path_to_data_dir, list_of_central_atoms)
print "total time is", time.time() - start_time
