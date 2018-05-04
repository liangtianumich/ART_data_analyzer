import matplotlib.pyplot as plt
import os
import pickle
import pandas as pd

displacement = pickle.load(open('disp_results.pkl','r'))
shear_strain = pickle.load(open('shear_strain_results.pkl','r'))
path_to_ovito = "ovito_min_sad_strain_disp"
ovito_strain = pd.read_csv(path_to_ovito, sep='\s+', skiprows = 1)
ovito_strain_xyz = ovito_strain.iloc[:,5:8]
ovito_strain_xyz.columns = ["vol_strain","shear_strain","disp"]

ovito_shear = []
ovito_vol = []
ovito_disp = []
my_shear = []
atom_item = []
my_disp = []

for index, row in ovito_strain_xyz.iterrows():
	ovito_shear.append(row["shear_strain"])
	ovito_vol.append(row["vol_strain"])
	ovito_disp.append(row["disp"])
	atom_item.append(index + 1)
	#row_xyz = ovito_strain_xyz.iloc[index]
	#ovito_xyz.append(row_xyz["shear_strain"])
	my_shear.append(shear_strain[index])
	my_disp.append(displacement[index])
plt.figure()
plt.plot(atom_item,ovito_shear,'ro',label="ovito_shear")
plt.plot(atom_item,my_shear,'bx',label="my shear")
plt.xlabel("atom item id")
plt.ylabel("shear strain")
plt.legend(loc='best')
plt.savefig("ovito_shear_strain_comparison.tif")
plt.close()

plt.figure()
plt.plot(atom_item,ovito_disp,'ro',label="ovito_disp")
plt.plot(atom_item,my_disp,'bx',label="my disp")
plt.xlabel("atom item id")
plt.ylabel("atomic displacement magnitude")
plt.legend(loc='best')
plt.savefig("ovito_displacement_comparison.tif")
plt.close()

	


