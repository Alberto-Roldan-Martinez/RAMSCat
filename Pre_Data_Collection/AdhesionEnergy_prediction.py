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
	popts = {													# popt predicted using Trend_AdhEnergy_Sites
                ('MgO', 'Au',  2, np.array([2.433770778924254, 1.0074054846113665, 2.18651222857168, 2.2455762462521727, 1.9020852732149833, 2.1746961810522736])),
                ('MgO', 'Au',  3, np.array([1.9542168619948455, 1.1014439587843892, 2.3084499403864833, 1.4544180750777675, 2.5022526376198715, 1.9706719475838992])),
                ('MgO', 'Au',  4, np.array([2.0141435446732086, 1.2811731386430516, 2.2745843866620783, 1.4222748682263509, 2.756306241456223, 1.9382905198186324])),
                ('MgO', 'Au',  5, np.array([1.9783290675048133, 1.3290659679379555, 2.2825255708618086, 1.4393564764211793, 3.2217670522828294, 1.817538002430754])),
				('MgO', 'Au',  6, np.array([2.319678310432216, 1.4245070028092464, 2.1961186288124823, 2.503923965931683, 3.13873938592232, 1.9138946081867072])),
				('MgO', 'Au',  7, np.array([2.132245747719973, 1.083640466904854, 2.2504750598392826, 1.5093307969890077, 2.897226198332409, 1.891606196681495])),
				('MgO', 'Au',  8, np.array([2.2600076123312363, 1.1253484003047012, 2.2319902426696445, 1.418855427818597, 3.3663694191918356, 1.886500500003469])),
				('MgO', 'Au',  9, np.array([2.167643459799678, 0.8860272081325079, 2.276217833042536, 1.6189887272725, 2.8145360818270233, 1.86491975732562])),
				('MgO', 'Au',  10, np.array([2.1438588251205704, 1.6363268933098492, 2.1586638828312874, 1.2692346875601228, 5.176499008423282, 1.5044917217215679])),
            }
	for i, sys in enumerate(popts):
		if sys[0] == support:
			if sys[1] == element:
				if sys[2] == icc:
					a, d_eq, r_eq, b, y_d_eq, y_r_eq = sys[3]
					i_adh_e = d_eq * (np.exp(-2*a*(x - r_eq)) - 2 * np.exp(-a*(x - r_eq*np.sin(y/x)))) +\
							  y_d_eq * (np.exp(-2*b*(y - y_r_eq)) - 2 * np.exp(-b*(y - y_r_eq*np.sin(y/x))))		# MORSE potentia
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
cluster_interface = {}
cluster_interface_cluster_neighbours = 0
for i in range(len(supported_cluster)):
	if supported_cluster[i].symbol in gas_cluster.get_chemical_symbols():
		cluster_indexes.append(supported_cluster[i].index)
	else:
		support_indexes.append(supported_cluster[i].index)
atom_cluster_neighbours = {}
for i in cluster_indexes:
	atom_cluster_neighbours[str(i)], atom_surface_neighbours = atom_neighbours(i, supported_cluster, cluster_indexes, support_indexes)
	if len(atom_surface_neighbours) > 0:
		cluster_interface_cluster_neighbours += len(atom_cluster_neighbours[str(i)])
		cluster_interface[str(i)] = atom_surface_neighbours

if len(cluster_interface) > 0:
	n_interface_cluster_atoms = len(cluster_interface)			# number of cluster atoms at the interface
#	print(cluster_interface, supported_cluster[i].symbol)
# list of unique coordinating atoms
#	unique_cluster_interface_indexes = [i for i in list(set(map(tuple, [cluster_interface[j] for j in cluster_interface])))[0]]
	average_cluster_coordination_interface_cluster_atoms = cluster_interface_cluster_neighbours/n_interface_cluster_atoms
else:
	n_interface_cluster_atoms = 0
	average_cluster_coordination_interface_cluster_atoms = 0


# Average of the shortest distances from atoms in the cluster to sites[0] in the support ==> along Z AXIS
average_shortest_cluster_site_distance = {}
z_c_interface = []
z_min_c_interface = min([supported_cluster.get_positions()[i][2] for i in cluster_indexes])
z_c_interface_indexes = [i for i in cluster_indexes if supported_cluster.get_positions()[i][2] < z_min_c_interface + 1]
z_max_surface = max([supported_cluster.get_positions()[i][2] for i in support_indexes])
z_s_neighbours_indexes = [i for i in support_indexes if supported_cluster.get_positions()[i][2] > z_max_surface - 1]
z_surface_neighbours = [supported_cluster.get_positions()[i][2] for i in z_s_neighbours_indexes]

for site in sites(support_name):
	for i in z_c_interface_indexes:
		z_c_interface.append(supported_cluster.get_positions()[i][2])
		distances = []
		shortest_cluster_site_distance = 0
		for j in z_s_neighbours_indexes:
			if supported_cluster[j].symbol == site:
				distances.append(supported_cluster.get_distance(i, j, mic=True, vector=False))
		shortest_cluster_site_distance += sorted(distances)[0]
	average_shortest_cluster_site_distance[str(site)] = shortest_cluster_site_distance / len(z_c_interface_indexes)

# Average distance between the support and the cluster interface atoms
average_z_interface = sum(z_c_interface) / len(z_c_interface)
average_z_surface = sum(z_surface_neighbours) / len(z_surface_neighbours)
average_interface_distance = average_z_interface - average_z_surface

# Predicting Adhesion energy
predicted_adhesion_e = []
for i in cluster_interface:
	isd_a = []
	isd_b = []
	for j in cluster_interface[str(i)]:
		if supported_cluster[j].symbol == sites(support_name)[0]:
			isd_a.append(supported_cluster.get_distance(i, j, mic=True, vector=False))
		elif supported_cluster[j].symbol == sites(support_name)[1]:
			isd_b.append(supported_cluster.get_distance(i, j, mic=True, vector=False))
	predicted_adhesion_e += morse_3D_Energies(support_name, supported_cluster[i].symbol,
											  len(atom_cluster_neighbours[str(i)]), sorted(isd_a)[0], sorted(isd_b)[0])

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
ifile = open("Trend_AdhEnergy.dat", 'w+')
ifile.write("#\t{}\n".format(support_name))
ifile.write("#\n# ic = n_interface_cluster_atoms\n# icc = average_cluster_coordination_interface_cluster_atoms\n")
ifile.write("# id = average distance from the cluster interface atoms to the surface neighbouring atoms\n")
ifile.write("# isd = average of the shortest distance from the interfacial atoms in the cluster to the surface site\n")
ifile.write("#\n# ic\ticc\tid\tisd_{}\tisd_{}\t\tE_Adh (eV)\tPredicted_E_Adh(eV)\tElements\tPath\n"
			.format(sites(support_name)[0], sites(support_name)[1]))
ifile.write("{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t{:>5.4f}\t\t{:>5.4f}\t{:>5.4f}" .format(n_interface_cluster_atoms,
																average_cluster_coordination_interface_cluster_atoms,
																		 average_interface_distance,
																		 average_shortest_cluster_site_distance[sites(support_name)[0]],
																		 average_shortest_cluster_site_distance[sites(support_name)[1]],
																		 adhesion_e, predicted_adhesion_e))
ifile.write("\t\t# {}\t\t{}\n" .format(list(set(gas_cluster.get_chemical_symbols()))[0], path_name))
ifile.close()
