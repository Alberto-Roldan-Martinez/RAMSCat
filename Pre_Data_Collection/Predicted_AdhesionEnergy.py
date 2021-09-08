'''
	Finds the cohesion energy of supported metallic systems

	USAGE: ~.py support_name(eg MgO) gas_cluster_path

'''

import os
import sys
import numpy as np
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


def morse_3D_Energies(support, element, icc, x, y):
	popts = {# support, metal, n-1											# popt and reference_e using Trend_AdhEnergy_Sites
#                ('MgO', 'Au',  1,  , 0.6396),	# 2 atoms ***
                ('MgO', 'Au',  2, 1.18296, 1.57988, 0.94593, 2.88265, 1.42953, 1.96775, 13.43070, 1.60617, -0.3442, -1.65847),		# 3 atoms
                ('MgO', 'Au',  3, 1.23430, 1.69203, 1.04552, 2.86776, 1.50758, 2.12566, 19.78013, 1.54962, -0.5365, -1.86308),		# 4 atoms
                ('MgO', 'Au',  4, 1.32852, 1.84366, 1.05141, 2.83982, 1.41161, 2.10295, 27.30569, 1.43322, -0.5228, -1.96177),		# 5 atoms
				('MgO', 'Au',  5, 1.36965, 1.80618, 1.18215, 2.65505, 1.70042, 2.12625, 12.25984, 1.63363,  0.2132, -1.91265),		# 6 atoms
				('MgO', 'Au',  6, 1.27389, 1.76288, 0.87434, 2.84338, 1.52271, 2.09858, 17.30035, 1.53972, -0.8052, -1.71807),		# 7 atoms
				('MgO', 'Au',  7, 1.25468, 1.77093, 0.86635, 2.89051, 1.45040, 2.10521, 24.11202, 1.44576, -0.8657, -1.90081),		# 8 atoms
				('MgO', 'Au',  8, 1.47385, 1.89203, 0.88521, 2.61610, 2.08310, 2.15869,  7.15395, 1.77023, -1.1498, -1.37068),		# 9 atoms
#				('MgO', 'Au',  9,  ),		# 10 atoms
				('MgO', 'Au', 10, 1.27326, 1.80154, 1.17582, 2.85178, 1.78804, 2.07628, 23.30049, 1.51652, -1.2182, -2.69864)		# 11 atom
            }
	for i, sys in enumerate(popts):
		if sys[0] == support:
			if sys[1] == element:
				if sys[2] == icc:
					support, element, icc, a1, a2, a_d_eq, a_r_eq, b1, b2, b_d_eq, b_r_eq, reference_e, e_min = sys
					i_adh_e = (a_d_eq * (np.exp(-2*a1*(x - a_r_eq)) - 2 * np.exp(-a2*(x - a_r_eq*np.sin(y/x)))) +
							   b_d_eq * (np.exp(-2*b1*(y - b_r_eq)) - 2 * np.exp(-b2*(y - b_r_eq*np.sin(y/x))))) + reference_e		# MORSE potentia

	return i_adh_e, reference_e, e_min


def atom_neighbours(atom_index, supported_cluster, cluster_indexes, support_indexes):
	cutoff = neighborlist.natural_cutoffs(supported_cluster, mult=1.3)
	a, b = neighborlist.neighbor_list('ij', supported_cluster, cutoff)
	atom_cluster_neighbours = sorted([b[n] for n in range(len(a)) if a[n] == atom_index and b[n] in cluster_indexes])
	atom_surface_neighbours = sorted([b[n] for n in range(len(a)) if a[n] == atom_index and b[n] in support_indexes])

	return atom_cluster_neighbours, atom_surface_neighbours


# Interface atoms and their coordination within the supported cluster
cluster_indexes = []
support_indexes = []
for i in range(len(supported_cluster)):
	if supported_cluster[i].symbol in gas_cluster.get_chemical_symbols():
		cluster_indexes.append(supported_cluster[i].index)
	else:
		support_indexes.append(supported_cluster[i].index)

# Atoms at the cluster and support inferfaces along Z AXIS
z_min_cluster_interface = min([supported_cluster.get_positions()[i][2] for i in cluster_indexes])
cluster_interface_indexes = [i for i in cluster_indexes if supported_cluster.get_positions()[i][2] <
							   z_min_cluster_interface + 1]
z_max_support = max([supported_cluster.get_positions()[i][2] for i in support_indexes])
support_neighbours_indexes = [i for i in support_indexes if supported_cluster.get_positions()[i][2] > z_max_support - 1]
support_neighbours = [supported_cluster.get_positions()[i][2] for i in support_neighbours_indexes]

# Average distance between the support and the cluster interface atoms
cluster_interface = []
for i in cluster_interface_indexes:
	cluster_interface.append([i, supported_cluster.get_positions()[i][2]])
cluster_interface.sort(key=lambda x: x[1])
if len(cluster_interface) > 1:
	average_z_interface = sum([i[1] for i in cluster_interface]) / len(cluster_interface)
else:
	average_z_interface = cluster_interface[0][1]
average_z_support = sum(support_neighbours) / len(support_neighbours)
average_interface_distance = average_z_interface - average_z_support

# Interface coordination
cluster_interface_adh_e = []
atom_cluster_neighbours = {}
atom_surface_neighbours = {}
average_shortest_cluster_site_distance = {}
average_shortest_cluster_site_distance[sites(support_name)[0]] = 0
average_shortest_cluster_site_distance[sites(support_name)[1]] = 0
supported_cluster.set_constraint()
cluster_interface_cluster_neighbours = 0
for i in cluster_interface_indexes:
	atom_cluster_neighbours[str(i)], atom_surface_neighbours[str(i)] = atom_neighbours(i, supported_cluster,
																					   cluster_indexes, support_indexes)
	distance_a = []
	distance_b = []
	shortest_cluster_site_distance = 0
	for j in support_neighbours_indexes:
		if supported_cluster[j].symbol == sites(support_name)[0]:
			distance_a.append([j, supported_cluster.get_distance(int(i), int(j), mic=True, vector=False)])
		elif supported_cluster[j].symbol == sites(support_name)[1]:
			distance_b.append([j, supported_cluster.get_distance(int(i), int(j), mic=True, vector=False)])
	distance_a.sort(key=lambda x: x[1])
	distance_b.sort(key=lambda x: x[1])
	i_atom_dz = supported_cluster.get_positions()[i][2] - average_z_support
	adh_e, reference_e, e_min = morse_3D_Energies(support_name, supported_cluster[int(i)].symbol, len(atom_cluster_neighbours[str(i)]),
									 distance_a[0][1], distance_b[0][1])

	cluster_interface_adh_e.append([i, round(adh_e, 5), round(e_min, 5), round(distance_a[0][1], 3), round(distance_b[0][1], 3)])

	if distance_a[0][0] not in atom_surface_neighbours[str(i)]:
		atom_surface_neighbours[str(i)].append(distance_a[0][0])
	if distance_b[0][0] not in atom_surface_neighbours[str(i)]:
		atom_surface_neighbours[str(i)].append(distance_b[0][0])
	average_shortest_cluster_site_distance[sites(support_name)[0]] += distance_a[0][1]/len(cluster_interface_indexes)
	average_shortest_cluster_site_distance[sites(support_name)[1]] += distance_b[0][1]/len(cluster_interface_indexes)
	cluster_interface_cluster_neighbours += len(atom_cluster_neighbours[str(i)])
average_cluster_coordination_interface_cluster_atoms = cluster_interface_cluster_neighbours/len(cluster_interface_indexes)
#unique_cluster_interface_indexes = [i for i in list(set(map(tuple, [cluster_interface[j] for j in cluster_interface])))[0]]

# Predict Adhesion energy
cluster_interface_adh_e.sort(key=lambda x: x[1])
primary_cluster_sites = []
secondary_cluster_sites = []
predicted_adhesion_e = 0
for n in range(len(cluster_interface_adh_e)):
	i = cluster_interface_adh_e[n][0]
	if i not in secondary_cluster_sites:
		primary_cluster_sites.append(i)
		for j in atom_cluster_neighbours[str(i)]:
			if j in cluster_interface_indexes:
				secondary_cluster_sites.append(j)

for n in range(len(cluster_interface_adh_e)):
# % margin of error in the minimum energy
	if cluster_interface_adh_e[0][1] < cluster_interface_adh_e[0][2]*0.85 and len(primary_cluster_sites) > 1:
		predicted_adhesion_e += cluster_interface_adh_e[n][1]/len(cluster_interface_adh_e)
		caca = 1
	elif cluster_interface_adh_e[0][1] < cluster_interface_adh_e[0][2]*0.75 and len(primary_cluster_sites) <= 1:
		if len(set(secondary_cluster_sites)) > 0:
			predicted_adhesion_e += cluster_interface_adh_e[n][1]/len(set(secondary_cluster_sites))
		else:
			predicted_adhesion_e += cluster_interface_adh_e[n][1]
		caca = 2
	elif cluster_interface_adh_e[n][0] in primary_cluster_sites:
		predicted_adhesion_e += cluster_interface_adh_e[n][1]/len(primary_cluster_sites)
		caca = 3
	elif cluster_interface_adh_e[n][0] in secondary_cluster_sites:
		predicted_adhesion_e += cluster_interface_adh_e[n][1]/len(set(secondary_cluster_sites))
		caca = 4

print(caca, predicted_adhesion_e, cluster_interface_adh_e, len(primary_cluster_sites), set(secondary_cluster_sites))#,"\n", [atom_surface_neighbours[str(i)] for i in cluster_interface_indexes], a_sites, len(a_sites))

# Calculate the array of neighbours coordination
neigh_coord = []
for i in cluster_interface_indexes:
	neigh_coord.append(len(atom_cluster_neighbours[str(i)]))

# Calculated Adhesion energy
clean_surf_file = str("~/RESEARCH/OTHER/DATASET/RPBE/Supports/"+support_name+"/"+support_name+"/Surface/"+input_file[0])
gas_cluster_file = str(gas_cluster_path + "/" + input_file[0])
supported_cluster_file = str("./" + input_file[0])
clean_surf = read(clean_surf_file)
gas_cluster = read(gas_cluster_file)
supported_cluster = read(supported_cluster_file)

e_supported = supported_cluster.get_total_energy()
e_surface = clean_surf.get_total_energy()
e_gas_cluster = gas_cluster.get_total_energy()
adhesion_e = e_supported - (e_gas_cluster + e_surface)

# Printing information
ifile = open("Predicted_AdhEnergy.dat", 'w+')
ifile.write("#\t{}\n".format(support_name))
ifile.write("#\n# ic1 = n primary interface_cluster_atoms\n# ic2 = n secondaty interface cluster atoms\n")
ifile.write("# icc = average_cluster_coordination_interface_cluster_atoms\n")
ifile.write("# id = average distance from the cluster interface atoms to the surface neighbouring atoms\n")
ifile.write("# isd = average of the shortest distance from the interfacial atoms in the cluster to the surface site\n")
ifile.write("#\n# ic1\t1c2\ticc\tid\tisd_{}\tisd_{}\t\tE_Adh (eV)\t\tCoordination [1..11]\t\t\tElements\tPath\n"
			.format(sites(support_name)[0], sites(support_name)[1]))
ifile.write("{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t\t{:>5.4f}\t{:>5.4f}\t  "
			.format(len(set(primary_cluster_sites)),	len(set(secondary_cluster_sites)),
					average_cluster_coordination_interface_cluster_atoms, average_interface_distance,
					average_shortest_cluster_site_distance[sites(support_name)[0]],
					average_shortest_cluster_site_distance[sites(support_name)[1]],	adhesion_e, predicted_adhesion_e))
for i in range(1, 12):
	ifile.write("{:>3d}".format(neigh_coord.count(i)))
ifile.write("\t\t# {}\t\t{}\n" .format(list(set(gas_cluster.get_chemical_symbols()))[0], path_name))
ifile.close()
