"""
	Versions:
		Alberto: 05/2022

"""

import os
from ase.io import read
from Coordination import Coordination
from GCN import Generalised_coodination
from Areas import Areas
from Zdistance import Cluster_surface_distance
from Energies import Energies
from WriteData import Write_labels, write_results


""" --------------------------- CLUSTER MODEL ---------------------------"""
cluster_elements = ["Au"]                           		# Elements in the Cluster
isolated_cluster = "./Cluster/OUTCAR"        				# file with the isolated cluster's energy
""" --------------------------- SURFACE MODEL---------------------------"""
support = "MgO"                             				# Support's name
support_size = [8, 8, 4]									# Dimension of the support's supercell
###################################################################################################

''' ---------------- Get and Print Results ---------------------'''
atoms = read("CONTCAR")
e_atoms = read("OUTCAR")
cluster = read(isolated_cluster)
coordination = Coordination(atoms, cluster_elements, support)
gcn = Generalised_coodination(atoms, cluster_elements, support)
area = Areas(atoms, cluster_elements, support)
z_distance = Cluster_surface_distance(atoms, cluster_elements, support)
energy = Energies(atoms, e_atoms, cluster_elements, cluster, support, support_size)

path = os.getcwd()
name = path.split("/")[-4]+"/"+path.split("/")[-3]+"/"+path.split("/")[-2]+"/"+path.split("/")[-1]

labels = ["N", "i_c", coordination.site_cluster_coordination_label, "i_cc", coordination.cluster_coord_labels,
				coordination.support_cluster_min_distance_labels, "cs_height", z_distance.zlabels, "GCN", "c_i_area",
				"c_s_area", "Esurf", "Ecoh", "Eadh", "Eb", "Etotal", "  structure_path"]
values = [coordination.cluster_size, coordination.interface_cluster, coordination.site_cluster_coordination,
		  coordination.interface_cc_average, coordination.cluster_ave_coordination, coordination.support_cluster_min_distance,
		  z_distance.interface_height, z_distance.cluster_cm_surface_distance, float(gcn.gcn_average),
		  area.cluster_interface_area, area.cluster_surface_area, energy.e_cluster_surface, energy.cohesion/coordination.cluster_size,
		  energy.adhesion, energy.binding/coordination.cluster_size, energy.e_total, name]

Write_labels("labels.txt", labels)
write_results("data.dat", values)


