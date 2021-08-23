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
                ('MgO', 'Au',  1, 1.8607086053315036, 1.208035121073158, 2.2261094470840157, 2.62965321327308, 1.576008693986123, 2.1578066234798876, 0.6396),	# 2 atoms ***
                ('MgO', 'Au',  2, 1.4927350400550732, 1.0624970024824945, 2.38020355343081, 1.4514764879189075, 2.282964315471724, 2.0120269550336984, 0.3442),		# 3 atoms
                ('MgO', 'Au',  3, 1.650710715231977, 1.2880965633594754, 2.2967050122715347, 1.4042094971557366, 2.214223474586348, 1.9386851604701019, 0.5365),	# 4 atoms
                ('MgO', 'Au',  4, 1.7042395478927133, 1.326251962144709, 2.3276014108681005, 1.6252060548342648, 3.471425276860701, 1.8848359939323371, 0.5228),	# 5 atoms
				('MgO', 'Au',  5, 2.0908097228375704, 1.4293216181936355, 2.2029959348956236, 2.6778975589780503, 3.230791785165967, 1.9099852159904098, 0.2051),	# 6 atoms
				('MgO', 'Au',  6, 1.7849198981431542, 1.0428222980300355, 2.271886577745904, 1.403426393718411, 2.448184836472611, 1.899406972802368, 0.8052),		# 7 atoms
				('MgO', 'Au',  7, 1.875042357550135, 1.004314385388471, 2.256700885018284, 1.407551007940099, 2.9138511768529005, 1.914316338871738, 0.8657),		# 8 atoms
				('MgO', 'Au',  8, 1.8599184127128634, 0.8864672354627985, 2.3061355646581925, 1.6044967310553975, 2.639307941186561, 1.8590276020279555, 1.1498),	# 9 atoms
				('MgO', 'Au',  9, 1.834262217337352, 1.742284321337311, 2.160124960457686, 1.1270893058451519, 4.710243723521435, 1.268619353701571, 1.0898),		# 10 atoms
				('MgO', 'Au',  10, )		# 11 atom
            }
	for i, sys in enumerate(popts):
		if sys[0] == support:
			if sys[1] == element:
				if sys[2] == icc:
					support, element, icc, a, d_eq, r_eq, b, y_d_eq, y_r_eq, reference_e = sys
					i_adh_e = d_eq * (np.exp(-2*a*(x - r_eq)) - 2 * np.exp(-a*(x - r_eq*np.sin(y/x)))) +\
							  y_d_eq * (np.exp(-2*b*(y - y_r_eq)) - 2 * np.exp(-b*(y - y_r_eq*np.sin(y/x)))) -\
							  reference_e		# MORSE potentia
	return i_adh_e


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
	cluster_interface.append(supported_cluster.get_positions()[i][2])
average_z_interface = sum(cluster_interface) / len(cluster_interface)
average_z_support = sum(support_neighbours) / len(support_neighbours)
average_interface_distance = average_z_interface - average_z_support

# Interface coordination
atom_cluster_neighbours = {}
average_shortest_cluster_site_distance = {}
average_shortest_cluster_site_distance[sites(support_name)[0]] = 0
average_shortest_cluster_site_distance[sites(support_name)[1]] = 0
supported_cluster.set_constraint()
cluster_interface_cluster_neighbours = 0
predicted_adhesion_e = 0
for i in cluster_interface_indexes:
	atom_cluster_neighbours[str(i)], atom_surface_neighbours = atom_neighbours(i, supported_cluster, cluster_indexes,
																			   support_indexes)
	distance_a = []
	distance_b = []
	shortest_cluster_site_distance = 0
	for j in support_neighbours_indexes:
		if supported_cluster[j].symbol == sites(support_name)[0]:
			distance_a.append(supported_cluster.get_distance(int(i), int(j), mic=True, vector=False))
		elif supported_cluster[j].symbol == sites(support_name)[1]:
			distance_b.append(supported_cluster.get_distance(int(i), int(j), mic=True, vector=False))
	average_shortest_cluster_site_distance[sites(support_name)[0]] += sorted(distance_a)[0]/len(cluster_interface_indexes)
	average_shortest_cluster_site_distance[sites(support_name)[1]] += sorted(distance_b)[0]/len(cluster_interface_indexes)
	predicted_adhesion_e += morse_3D_Energies(support_name, supported_cluster[int(i)].symbol,
									len(atom_cluster_neighbours[str(i)]), sorted(distance_a)[0], sorted(distance_b)[0])
	cluster_interface_cluster_neighbours += len(atom_cluster_neighbours[str(i)])

#predicted_adhesion_e = predicted_adhesion_e / len(cluster_interface_indexes)


average_cluster_coordination_interface_cluster_atoms = cluster_interface_cluster_neighbours/len(cluster_interface_indexes)
#unique_cluster_interface_indexes = [i for i in list(set(map(tuple, [cluster_interface[j] for j in cluster_interface])))[0]]

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
adhesion_e = (e_supported - (e_gas_cluster + e_surface))

# Printing information
ifile = open("Predicted_AdhEnergy.dat", 'w+')
ifile.write("#\t{}\n".format(support_name))
ifile.write("#\n# ic = n_interface_cluster_atoms\n# icc = average_cluster_coordination_interface_cluster_atoms\n")
ifile.write("# id = average distance from the cluster interface atoms to the surface neighbouring atoms\n")
ifile.write("# isd = average of the shortest distance from the interfacial atoms in the cluster to the surface site\n")
ifile.write("#\n# ic\ticc\tid\tisd_{}\tisd_{}\t\tE_Adh (eV)\t\tElements\tPath\n"
			.format(sites(support_name)[0], sites(support_name)[1]))
ifile.write("{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t\t{:>5.4f}\t{:>5.4f}" .format(len(cluster_interface_indexes),
																average_cluster_coordination_interface_cluster_atoms,
																		 average_interface_distance,
																		 average_shortest_cluster_site_distance[sites(support_name)[0]],
																		 average_shortest_cluster_site_distance[sites(support_name)[1]],
																		 adhesion_e, predicted_adhesion_e))
ifile.write("\t\t# {}\t\t{}\n" .format(list(set(gas_cluster.get_chemical_symbols()))[0], path_name))
ifile.close()
