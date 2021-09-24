'''
	Versions:
		Alberto: 08/2019

	STRUCTURE:
		- check structure and energy (OUTCAR) files
		- check boundary conditions
		- Import coordination data
		- Import energy from the OUTCAR

'''

import os, sys
from CheckFiles import checkFiles
from Coordination import Coordination
from GCN import Generalised_coodination
from Areas import Areas
from Zdistance import Cluster_surface_distance
from Energies import Energies
from WriteData import Write_labels, write_results

if sys.argv:
	argument = sys.argv[1]
else:
	argument = "."

path = os.getcwd()
name = path.split("/")[-4]+"/"+path.split("/")[-3]+"/"+path.split("/")[-2]+"/"+path.split("/")[-1]

cluster_elements = ["Au"]                           		# Elements in the Cluster
support = "MgO"                             				# Support's name
inputfiles = ["OUTCAR", "CONTCAR"]							# files containing the energies and geometries
isolated_cluster = argument + "/gas/OUTCAR"           	# file with the isolated cluster's energy
isolated_support = "/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Supports/MgO/MgO/Surface/OUTCAR"	# file with the support's energy


coordination = Coordination(inputfiles[1], cluster_elements, support)
gcn = Generalised_coodination(inputfiles[1], cluster_elements, support)
area = Areas(inputfiles[1], cluster_elements, support)
z_distance = Cluster_surface_distance(inputfiles[1], cluster_elements, support)
energy = Energies(inputfiles, isolated_support, isolated_cluster, cluster_elements, support)


labels = ["N", "i_c", coordination.site_cluster_coordination_label, "i_cc", coordination.cluster_coord_labels,
				coordination.support_cluster_min_distance_labels, "cs_height", z_distance.zlabels, "GCN", "c_i_area",
				"c_s_area", "Esurf", "Ecoh", "Eadh", "Eb", "Etotal", "  structure_path"]
values = [coordination.cluster_size, coordination.interface_cluster, coordination.site_cluster_coordination,
		  coordination.interface_cc_average, coordination.cluster_ave_coordination, coordination.support_cluster_min_distance,
		  z_distance.interface_height, z_distance.cluster_cm_surface_distance, float(gcn.gcn_average),
		  area.cluster_interface_area, area.cluster_surface_area, energy.e_cluster_surface, energy.cohesion,
		  energy.adhesion, energy.binding, energy.e_total, name]

Write_labels("labels.txt", labels)
write_results("data.dat", values)


