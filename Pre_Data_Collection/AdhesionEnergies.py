'''
	Finds the cohesion energy of supported metallic systems

	USAGE: ~.py support_name(eg MgO) gas_cluster_path

'''

import os
import sys
from ase.io import read
from ase import neighborlist
from Library import sites


input_file = ["OUTCAR", "CONTCAR"]
support_name = sys.argv[1] 			# name of clean support, e.g., MgO
gas_cluster_path = sys.argv[2]		# path to the isolated (gas phase) cluster with exact same geometry to the supported
path = os.getcwd()
path_name = path.split("/")[-4]+"/"+path.split("/")[-3]+"/"+path.split("/")[-2]+"/"+path.split("/")[-1]

#        CHECK for information required
try:
	clean_surf_file = str("~/RESEARCH/OTHER/DATASET/RPBE/Supports/"+support_name+"/"+support_name+"/Surface/"+input_file[1])
	clean_surf = read(clean_surf_file)
except:
	print(" Clean Support Name Incorrect or Inexistent!")
	exit()
try:
	gas_cluster = read(gas_cluster_path)
except:
	gas_cluster_file = str(gas_cluster_path + "/" + input_file[1])
	gas_cluster = read(gas_cluster_file)
try:
	supported_cluster_file = str("./" + input_file[1])
	supported_cluster = read(supported_cluster_file, index=-1)
except:
	print(" Supported Cluster Input does not exist in the current directory!")


def atom_neighbours(atom_index, supported_cluster, cluster_indexes, support_indexes):
	cutoff = neighborlist.natural_cutoffs(supported_cluster, mult=1.15)
	a, b = neighborlist.neighbor_list('ij', supported_cluster, cutoff)
	atom_cluster_neighbours = [b[n] for n in range(len(a)) if a[n] == atom_index and b[n] in cluster_indexes]
	atom_surface_neighbours = [b[n] for n in range(len(a)) if a[n] == atom_index and b[n] in support_indexes]
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
	atom_cluster_neighbours, atom_surface_neighbours = atom_neighbours(i, supported_cluster, cluster_indexes, support_indexes)
	if len(atom_surface_neighbours) > 0:
		cluster_interface_cluster_neighbours += len(atom_cluster_neighbours)
		cluster_interface[str(i)] = atom_surface_neighbours

if len(cluster_interface) > 0:
	n_interface_cluster_atoms = len(cluster_interface)			# number of cluster atoms at the interface
	average_cluster_coordination_interface_cluster_atoms = cluster_interface_cluster_neighbours/n_interface_cluster_atoms
else:
	n_interface_cluster_atoms = 0
	average_cluster_coordination_interface_cluster_atoms = 0


# Average of the shortest distances from atoms in the cluster to sites[0] in the support ==> along Z AXIS
average_shortest_cluster_site_distance = 0
z_c_interface = []
z_min_c_interface = min([supported_cluster.get_positions()[i][2] for i in cluster_indexes])
z_c_interface_indexes = [i for i in cluster_indexes if supported_cluster.get_positions()[i][2] < z_min_c_interface + 1]
z_max_surface = max([supported_cluster.get_positions()[i][2] for i in support_indexes if supported_cluster[i].symbol ==
					 sites(support_name)[0]])
z_s_neighbours_indexes = [i for i in support_indexes if supported_cluster.get_positions()[i][2] > z_max_surface - 1
							  and supported_cluster[i].symbol == sites(support_name)[0]]
z_surface_neighbours = [supported_cluster.get_positions()[i][2] for i in z_s_neighbours_indexes]

for i in z_c_interface_indexes:
	z_c_interface.append(supported_cluster.get_positions()[i][2])
	distances = []
	for j in z_s_neighbours_indexes:
		distances.append([supported_cluster[j].symbol, supported_cluster.get_distance(i, j, mic=True, vector=False)])
	distances.sort(key=lambda x: x[1])
	average_shortest_cluster_site_distance += distances[0][1]
average_shortest_cluster_site_distance = average_shortest_cluster_site_distance / len(z_c_interface_indexes)

# Average distance between the support and the cluster interface atoms
average_z_interface = sum(z_c_interface) / len(z_c_interface)
average_z_surface = sum(z_surface_neighbours) / len(z_surface_neighbours)
average_interface_distance = average_z_interface - average_z_surface

# Adhesion energy
clean_surf_file = str("~/RESEARCH/OTHER/DATASET/RPBE/Supports/"+support_name+"/"+support_name+"/Surface/"+input_file[0])
gas_cluster_file = str(gas_cluster_path + "/" + input_file[0])
supported_cluster_file = str("./" + input_file[0])
clean_surf = read(clean_surf_file)
gas_cluster = read(gas_cluster_file)
supported_cluster = read(supported_cluster_file)

e_supported = supported_cluster.get_total_energy()
e_surface = clean_surf.get_total_energy()
e_gas_cluster = gas_cluster.get_total_energy()
adhesion_e = (e_supported - (e_gas_cluster + e_surface))

# Printing information
ifile = open("Trend_AdhEnergy.dat", 'w+')
ifile.write("#\n# ic = n_interface_cluster_atoms\n# icc = average_cluster_coordination_interface_cluster_atoms\n")
ifile.write("# id = average distance from the cluster interface atoms to the surface neighbouring atoms\n")
ifile.write("# isd = average of the shortest distance from the interfacial atoms in the cluster to the surface most favourable site (in Library)\n")
ifile.write("#\n# ic\ticc\tid\tisd\t\tE_Adh (eV)\tElements\tPath\n")
ifile.write("{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t\t{:>5.4f}" .format(n_interface_cluster_atoms,
																average_cluster_coordination_interface_cluster_atoms,
																		 average_interface_distance,
																		 average_shortest_cluster_site_distance,
																		 adhesion_e))
ifile.write("\t\t# {}\t\t{}\n" .format(list(set(gas_cluster.get_chemical_symbols()))[0], path_name))
ifile.close()
