#!/usr/bin/env python
"""
this module runs ART simulations in slurm cluster environment,
the python multiprocessing module only works on single node so that
we will submit the job to compute node with at least 64 cores
"""
import subprocess
import os
import shutil
def run_art_cluster_slurm(path_to_data_dir, input_param):
	"""
	this function will invoke running artn simulations in cluster set up by 
	submitting a job to one compute node, which performs a parallelized
	art running using multiprocessing module with its maximum cores.
	
	Note:
		the key "num_of_proc" in input setting file will be interpreted as 
	the number of compute nodes when in this cluster mode
	"""
	#path_to_input_files = input_param['path_to_input_files']
	num_of_proc = input_param["num_of_proc"]
	central_atom_list = input_param['central_atom_list']
	
	# copy these two files into another dir in path_to_data_dir
	path_to_submit_slurm = os.path.join(os.environ["MY_ART"], "src/python/ART_wrapper/submit_slurm.sh")
	path_to_job_slurm = os.path.join(os.environ["MY_ART"], "src/python/ART_wrapper/job_slurm.py")
	path_to_slurm_output = os.path.join(path_to_data_dir, "slurm")
	if not os.path.isdir(path_to_slurm_output):
		os.makedirs(path_to_slurm_output)
	
	for full_file_name in [path_to_submit_slurm, path_to_job_slurm]:
		if os.path.isfile(full_file_name):
			shutil.copy(full_file_name,path_to_slurm_output)
	
	path_to_sh_file = os.path.join(path_to_slurm_output,"submit_slurm.sh")
	path_to_job_file = os.path.join(path_to_slurm_output,"job_slurm.py")
	# split the central_atom_list into the num_of_proc folds
	all_tests_folds = list(split(central_atom_list, num_of_proc))
	for test_fold in all_tests_folds:
		print "current job central atoms list:", test_fold
		test_fold = str(test_fold).replace(",", "")
		test_fold = "'%s'"%test_fold		
		subprocess.check_call("export ATOM_LIST=%s;export PATH_TO_JOB=%s;sbatch %s"%(test_fold,path_to_job_file,path_to_sh_file), shell=True)
		#run_bart_node(list_of_run_dir)
	print "done submitting jobs of running ART simulations using %s of compute nodes!"%num_of_proc
	
def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))
