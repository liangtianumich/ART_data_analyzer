import zipfile
import os
import json
from util import data_dir_to_test_dir,update_central_atom_list

def archive_project(path_to_data_dir):
	"""
	this function archive the most necessary data inside each test of an art data project
	for researchers to share data between each other
	
	First, it will archive the final_selected_events.json, which specifies the
	related events data to be saved
	
	Second, it will archive all configuration files that are related to the final filtered events
	saved in the final_selected_events.json
	
	Finally, it will also archive the log.file.1, events.list, bart.sh or mod_bart.sh file
	that are needed for many other functions
	"""
	print 'creating archive art_data_project.zip in %s'%path_to_data_dir
	path_to_zip = os.path.join(path_to_data_dir,'art_data_project.zip')
	# compress the file size
	zf = zipfile.ZipFile(path_to_zip,'w',zipfile.ZIP_DEFLATED,allowZip64 = True)
	
	print "changing directory to %s"%path_to_data_dir
	os.chdir(path_to_data_dir)
	
	list_of_files = archived_file_names(path_to_data_dir)
	#operation_on_events(path_to_data_dir, list_of_test_id, operation, num_of_proc)
	#for file_name in list_of_file_names:
	#operation = lambda x: add_file_name(x,zf)
	#result_list = pool.map(partial(add_file_name,zf = zf),list_of_files)
	for file_name in list_of_files:
		add_file_name(file_name,zf)

def archived_file_names(path_to_data_dir):
	path_to_final_selected_events = os.path.join(path_to_data_dir,"final_selected_events.json")
	#path_to_central_atom_list = os.path.join(path_to_data_dir,"central_atom_list.json")
	#path_to_int_atom_list = os.path.join(path_to_data_dir, "interested_atom_list.json")
	
	print "reading final_selected_events.json"
	if os.path.exists(path_to_final_selected_events):
		print "loading final_selected_events.json, ensure that you always get most updated final_selected_events.json by --filter --re_calc if you have calculated more tests"
		final_selected_events = json.load(open(path_to_final_selected_events,"r"))
	else:
		raise Exception("final_selected_events.json does not exist in %s"%path_to_data_dir)

	list_of_file_names = [["final_selected_events.json",path_to_final_selected_events]]
	
	#if os.path.exists(path_to_central_atom_list):
	#	print "archiving central_atom_list.json, ensure central_atom_list.json is updated by --update_input if using --art --run_more"
	#	central_file = ["central_atom_list.json",path_to_central_atom_list]
	#	list_of_file_names.append(central_file)
	
	#if os.path.exists(path_to_int_atom_list):
	#	print "archiving interested_atom_list.json"
	#	int_file = ["interested_atom_list.json",path_to_int_atom_list]
	#	list_of_file_names.append(int_file)
	file_dirs = os.listdir(path_to_data_dir)
	for f in file_dirs:
		if f == "final_selected_events.json" or f == "art_data_project.zip" or f == "central_atom_list.json":
			continue
		path_to_file = os.path.join(path_to_data_dir,f)
		if os.path.isfile(path_to_file):
			list_of_file_names.append([f,path_to_file])
			
	all_tests_id = []
	for event in final_selected_events:
		test_id = int(event[0][4:])
		all_tests_id.append(test_id)
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
	
	final_tests_id = list(set(all_tests_id))
	# update the central_atom_list.json
	print "update the central_atom_list.json to be all test id saved in final_selected_events.json"
	update_central_atom_list(path_to_data_dir,final_tests_id)
	list_of_file_names.append(["central_atom_list.json", os.path.join(path_to_data_dir, "central_atom_list.json")])
	
	for test in final_tests_id:
		path_to_test_dir = data_dir_to_test_dir(path_to_data_dir,test)
		test_file = return_test_results_file_name(path_to_test_dir)
		list_of_file_names.extend(test_file)
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

def return_test_results_file_name(path_to_test_dir):
	list_of_results_files = ["log.file.1","events.list","bart.sh","mod_bart.sh"]
	file_name_path = []
	for curr_file in list_of_results_files:
		path_to_file = os.path.join(path_to_test_dir,curr_file)
		if os.path.isfile(path_to_file):
			file_name = "/%s/%s"%(os.path.basename(path_to_test_dir),curr_file)
			file_name_path.append([file_name,path_to_file])
		#raise Exception("in %s, file %s does not exist"%(path_to_test_dir,curr_file))
	return file_name_path
	
def add_file_name(file_names,zf):
	file_name,path_to_file = file_names[0],file_names[1]
	print 'adding file %s'%path_to_file
	zf.write(path_to_file,file_name)
