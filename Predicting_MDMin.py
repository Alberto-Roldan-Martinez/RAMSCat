#!/usr/bin/env python3
#"""
#    Versions:
#        Alberto: 08/2019
#        Alberto: 09/2020
#
#"""

import os, sys
from ase.io import read
import ase.io.vasp
from ase.optimize import MDMin
from RAMSCat import RAMSCat
from Coordination import Coordination, Generalised_coodination
from Energies import Energy_prediction
from Properties import Areas, Cluster_surface_distance
from WriteData import Write_labels, write_results, write_out

#####################################################################################################
cluster_elements = [i for i in sys.argv[1].split("-")]                      # Elements in the Cluster
structurefile = sys.argv[2] 												# file name, e.g. POSCAR
support = sys.argv[3]                                                       # Surface name
support_size = [sys.argv[4], sys.argv[5], sys.argv[6]] #"/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Supports/MgO/MgO/Surface/OUTCAR"
fmax = float(sys.argv[7])
####################################################################################################

path = os.getcwd()
name = path.split("/")[-4]+"/"+path.split("/")[-3]+"/"+path.split("/")[-2]+"/"+path.split("/")[-1]

''' --------------- Structure Optimisation ---------------------'''
atoms = read(structurefile)
atoms.calc = RAMSCat(atoms, cluster_elements, support, support_size)

dyn = MDMin(atoms, logfile='Optimisation.txt', trajectory='trajectory.traj')
dyn.run(fmax=fmax) #5, steps=150)
ase.io.vasp.write_vasp("Optimised.vasp", atoms, direct=False, vasp5=True, sort=True, ignore_constraints=False)

''' ---------------- Get and Print Results ---------------------'''
coordination = Coordination(atoms, cluster_elements, support)
gcn = Generalised_coodination(atoms, cluster_elements, support)
area = Areas(atoms, cluster_elements, support)
z_distance = Cluster_surface_distance(atoms, cluster_elements, support)
energies = Energy_prediction(atoms, cluster_elements, support, support_size)

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
