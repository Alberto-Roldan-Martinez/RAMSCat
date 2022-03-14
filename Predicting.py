#!/usr/bin/env python3
#"""
#    Versions:
#        Alberto: 08/2019
#        Alberto: 09/2020
#
#"""

import os, sys
from Coordination import Coordination
from GCN import Generalised_coodination
from Areas import Areas
from Zdistance import Cluster_surface_distance
from Energies import Energy_prediction
from WriteData import Write_labels, write_results, write_out

#####################################################################################################
cluster_elements = [i for i in sys.argv[1].split("-")]                      # Elements in the Cluster
structurefile = sys.argv[2] 												# file name, e.g. POSCAR
support = sys.argv[3]                                                       # Surface name
support_size = sys.argv[4] #"/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Supports/MgO/MgO/Surface/OUTCAR"
####################################################################################################

path = os.getcwd()
name = path.split("/")[-4]+"/"+path.split("/")[-3]+"/"+path.split("/")[-2]+"/"+path.split("/")[-1]

coordination = Coordination(structurefile, cluster_elements, support)
gcn = Generalised_coodination(structurefile, cluster_elements, support)
area = Areas(structurefile, cluster_elements, support)
z_distance = Cluster_surface_distance(structurefile, cluster_elements, support)
energies = Energy_prediction(structurefile, cluster_elements, support, support_size)

labels = ["N", "i_c", coordination.site_cluster_coordination_label, "i_cc", coordination.cluster_coord_labels,
				coordination.support_cluster_min_distance_labels, "cs_height", z_distance.zlabels, "GCN", "c_i_area",
				"c_s_area", "Esurf", "Ecoh", "Eadh", "Eb", "Etotal", "  structure_path"]
values = [coordination.cluster_size, coordination.interface_cluster, coordination.site_cluster_coordination,
		  coordination.interface_cc_average, coordination.cluster_ave_coordination, coordination.support_cluster_min_distance,
		  z_distance.interface_height, z_distance.cluster_cm_surface_distance, float(gcn.gcn_average),
		  area.cluster_interface_area, area.cluster_surface_area, energies.e_cluster_surface, energies.e_coh/coordination.cluster_size,
		  energies.e_adh, energies.e_binding/coordination.cluster_size, energies.e_total, name]


Write_labels("Predicted.txt", labels)
write_results("Predicted.dat", values)
os.system("cat Predicted.dat >> Predicted.txt; mv Predicted.txt Predicted.dat")
write_out(structurefile, energies.e_total)

