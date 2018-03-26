This ART data analyzer package is a python package that analyze the initial and saddle point configurations of each possible event, output the information such as strain (e.g. the volumetric and shear strain) as a statistical distribution etc. This data package can be decomposed and organized into the following modules:

Data modules: extracting data from various data file formats (default format without extension) , output data into 
initial configurations, saddle point configurations, final configurations, with option to not extract final state configurations, with configurations as a class object. 
Class Configuration(..):
 def __init__(atom_num,position_df):
   Self.atom_num = atom_num
   Self. position_df = pd.DataFrame( index, column=)

 def data_to_df():





 

Current version 1.00 support:
Strain analyze, lammps data

Operating system:
This python package tested on Mac.OS.Sierra.10.12.6

How to install this python package:
This package is written in python scripting language to analyze the ART data in various file format that user can specify; there is no need to compile and build.  You can put this package in any dir. But you need to source the environment.sh file before using this package, the purpose is to create the necessary environmental variables for current terminal/bash sessions for the python interpreter/bash to find the probable PYTHONPATH/PATH. Some environmental variables are also needed to point to the specific dir path for finding the correct python modules.

How to use this python package:
/src contains all the source code of various python libraries, it may also contains other code such as Matlab code for the python to run using subprocess module, these code will be store inside /src/matlab

/scripts: contains executable scripts that you will directly invoke on the command line
/examples: contains examples of executable running scripts and demo data. 

Executable files:
art_data:
Prompt the options that you can use for various analyze feature of the python script
Use sys.parser to parse the user input

-calc/-c invoke the art_data_calculator exe 
-v invoke the visualization exe
-m invoke the correlation model exe by adopting the various data models implemented in python scikit-learn machine learning package
…





Art_data_


 


Features/modules available in this package

event_locator(): locate the event to be analyzed


strain_calculator_event(): calculate the strain for this specific event


strain calculator(): calculate various strain events

