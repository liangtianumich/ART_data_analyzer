import zipfile
import os
import json
import multiprocessing as mp
from functools import partial
from util import data_dir_to_test_dir

def archive_project(path_to_data_dir):
	
	print 'creating archive art_data_project.zip in %s'%path_to_data_dir
	path_to_zip = os.path.join(path_to_data_dir,'art_data_project.zip')
	zf = zipfile.ZipFile(path_to_zip, mode='w')
	
	print "changing directory to %s"%path_to_data_dir
	os.chdir(path_to_data_dir)
	
	#num_of_proc = input_param["num_of_proc"]
	
	#pool = mp.Pool(processes = num_of_proc)
	
	list_of_files = archived_file_names(path_to_data_dir)
	#operation_on_events(path_to_data_dir, list_of_test_id, operation, num_of_proc)
	#for file_name in list_of_file_names:
	#operation = lambda x: add_file_name(x,zf)
	#result_list = pool.map(partial(add_file_name,zf = zf),list_of_files)
	for file_name in list_of_files:
		add_file_name(file_name,zf)

def archived_file_names(path_to_data_dir):
	path_to_final_selected_events = os.path.join(path_to_data_dir,"final_selected_events.json")
	print "reading final_selected_events.json"
	final_selected_events = json.load(open(path_to_final_selected_events,"r"))
	list_of_file_names = [["final_selected_events.json",path_to_final_selected_events]]
	for event in final_selected_events:
		test_id = int(event[0][4:])
		path_to_test_dir = data_dir_to_test_dir(path_to_data_dir,test_id)
		init, sad, fin = event[1][0],event[1][1],event[1][2]
		#path_to_init = os.path.join(path_to_test_dir,init)
		#path_to_sad = os.path.join(path_to_test_dir,sad)
		#path_to_fin = os.path.join(path_to_test_dir,fin)
		init_file = return_file_names(path_to_test_dir, init)
		sad_file = return_file_names(path_to_test_dir, sad)
		fin_file = return_file_names(path_to_test_dir, fin)
		list_of_file_names.append(init_file)
		list_of_file_names.append(sad_file)
		list_of_file_names.append(fin_file)
	return list_of_file_names		
		
def return_file_names(path_to_test_dir, init):
	path_to_init = os.path.join(path_to_test_dir,init)
	if os.path.isfile(path_to_init):
		file_name = "/%s/%s"%(os.path.basename(path_to_test_dir),init)
		return [file_name,path_to_init]
	elif os.path.isfile(path_to_init+".dump"):
		file_name = "/%s/%s"%(os.path.basename(path_to_test_dir),init+".dump")
		return [file_name,path_to_init+".dump"]
	else:
		raise Exception("in %s, has no file named %s or %s.dump"%(path_to_test_dir,init,init))	
	
def add_file_name(file_names,zf):
	file_name,path_to_file = file_names[0],file_names[1]
	print 'adding file %s'%path_to_file
	zf.write(path_to_file,file_name)
