'''
	Finds the cohesion energy of supported metallic systems

	USAGE: ~.py path_for_the_DFT_output_files_contains_the_energy_and_geometry

'''

import os
import sys
from ase.io import read
from ase import neighborlist
from Library import isolated_atoms, ecoh_bulk


input_file = "OUTCAR"
support_name = sys.argv[1] 			# name of clean support, e.g., MgO
gas_cluster_path = sys.argv[2]		# path to the isolated (gas phase) cluster with exact same geometry to the supported
path = os.getcwd()
path_name = path.split("/")[-4]+"/"+path.split("/")[-3]+"/"+path.split("/")[-2]+"/"+path.split("/")[-1]

#        CHECK for information required
try:
	clean_surf_file = str("~/RESEARCH/OTHER/DATASET/RPBE/Supports/"+support_name+"/"+support_name+"/Surface/"+input_file)
	clean_surf = read(clean_surf_file)
except:
	print(" Clean Support Name Incorrect or Inexistent!")
	exit()
try:
	gas_cluster = read(gas_cluster_path)
except:
	gas_cluster_file = str(gas_cluster_path + "/" + input_file)
	gas_cluster = read(gas_cluster_file)
try:
	supported_cluster_file = str("./" + input_file)
	supported_cluster = read(supported_cluster_file)
except:
	print(" Supported Cluster Input does not exist in the current directory!")

def atom_neighbours(atom_index, supported_cluster, cluster_indexes):
	cutoff = neighborlist.natural_cutoffs(supported_cluster, mult=1.25)
	a, b = neighborlist.neighbor_list('ij', supported_cluster, cutoff)
	atom_cluster_neighbours = [b[i] for i in range(len(a)) if a[i] == atom_index and b[i] in cluster_indexes]
	atom_surface_neighbours = [b[i] for i in range(len(a)) if a[i] == atom_index and b[i] not in cluster_indexes]
	return atom_cluster_neighbours, atom_surface_neighbours

# Interface atoms and their coordination within the supported cluster
cluster_indexes = []
support_indexes = []
cluster_interface = {}
cluster_interface_cluster_neighbours = 0
for i in range(len(supported_cluster)):
	if supported_cluster[i].symbol in gas_cluster.get_chemical_symbols():
		cluster_indexes.append(supported_cluster[i].index)
	else:
		support_indexes.append(supported_cluster[i].index)
for i in cluster_indexes:
	atom_cluster_neighbours, atom_surface_neighbours = atom_neighbours(i, supported_cluster, cluster_indexes)
	if len(atom_surface_neighbours) > 0:
		cluster_interface_cluster_neighbours += len(atom_cluster_neighbours)
		cluster_interface[str(i)] = atom_surface_neighbours

n_interface_cluster_atoms = len(cluster_interface)			# number of cluster atoms at the interface
average_cluster_coordination_interface_cluster_atoms = cluster_interface_cluster_neighbours/n_interface_cluster_atoms


# CORRET


# Cluster_interface distance to the surface
interface_distance = 0
for i in cluster_interface:
	interface_distance += sorted([supported_cluster.get_distance(i, j, mic=True, vector=False) for j in cluster_interface[str(i)]])[0]
average_interface_distance = interface_distance / len(cluster_interface)



# Adhesion energy
e_supported = supported_cluster.get_total_energy()
e_surface = clean_surf.get_total_energy()
e_gas_cluster = gas_cluster.get_total_energy()

adhesion_e = (e_supported - (e_gas_cluster + e_surface))

# Printing information
ifile = open("Trend_AdhEnergy.dat", 'w+')
ifile.write("# ic = n_interface_cluster_atoms\n# icc = average_cluster_coordination_interface_cluster_atoms\n")
ifile.write("# id = average distance from the cluster interface atoms to the closest surface atoms\n")
ifile.write("# ic\ticc\tid\tE_Adh (eV)\t\tElements\tPath\n")
ifile.write("{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t\t{:>5.4f}\t" .format(n_interface_cluster_atoms,
														   average_cluster_coordination_interface_cluster_atoms,
														   average_interface_distance, adhesion_e))
ifile.write("\t# {}\t\t{}\n" .format(set(gas_cluster.get_chemical_symbols()), path_name))
ifile.close()
