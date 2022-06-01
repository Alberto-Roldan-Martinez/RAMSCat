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
from ase.optimize import BFGS, MDMin
from RAMSCat import RAMSCat
from Coordination import Coordination, Generalised_coodination
from Energies import Energy_prediction
from Properties import Areas, Cluster_surface_distance, Mean_interatomic_distance, Sphericity
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

dyn = BFGS(atoms, logfile='Optimisation.txt') #, trajectory='trajectory.traj')
##dyn = MDMin(atoms, logfile='Optimisation.txt') #, trajectory='trajectory.traj')
#dyn.run(fmax=fmax) #, steps=500)
ase.io.vasp.write_vasp("CONTCAR.vasp", atoms, direct=False, vasp5=True, sort=True, ignore_constraints=False)

''' ---------------- Get and Print Results ---------------------'''
values = [Coordination(atoms, cluster_elements, support).cluster_size,  							# N
          Coordination(atoms, cluster_elements, support).interface_cluster,  						# i_c
          Coordination(atoms, cluster_elements, support).site_cluster_coordination,  				# site(s)
          Coordination(atoms, cluster_elements, support).interface_cc_average,  					# i_cc
          Coordination(atoms, cluster_elements, support).cluster_ave_coordination,  				# cc
          float(Generalised_coodination(atoms, cluster_elements, support).gcn_average),  			# GCN
          Coordination(atoms, cluster_elements, support).support_cluster_min_distance,  			# dist_X
          Cluster_surface_distance(atoms, cluster_elements, support).interface_height,  			# cs_dist
          Cluster_surface_distance(atoms, cluster_elements, support).cluster_cm_surface_distance,  	# cm_dist
          Mean_interatomic_distance(atoms, cluster_elements, support).mean_distance,  				# cc_dist
          Sphericity(atoms, cluster_elements, support).shape_ratio,  								# shape
          Sphericity(atoms, cluster_elements, support).sphericity,  								# sphericity
          Areas(atoms, cluster_elements, support).cluster_interface_area,  							# c_i_area
          Areas(atoms, cluster_elements, support).cluster_surface_area,  							# c_s_area
          Energy_prediction(atoms, cluster_elements, support, support_size).e_cluster_surface,
          Energy_prediction(atoms, cluster_elements, support, support_size).cohesion,
          Energy_prediction(atoms, cluster_elements, support, support_size).adhesion,
          Energy_prediction(atoms, cluster_elements, support, support_size).binding,
          Energy_prediction(atoms, cluster_elements, support, support_size).e_total,
          name]

labels = ["N", "i_c", Coordination(atoms, cluster_elements, support).site_cluster_coordination_label, "i_cc",
          Coordination(atoms, cluster_elements, support).cluster_coord_labels, "GCN",
          Coordination(atoms, cluster_elements, support).support_cluster_min_distance_labels, "cs_dist", "cm_dist",
          "cc_dist", "shape", "sphericity", "c_i_area", "c_s_area", "Esurf", "Ecoh", "Eadh", "Eb", "Etotal",
          "  structure_path"]

Write_labels("Predicted.tmp", labels)
write_results("Predicted.dat", labels, values)
os.system("cat Predicted.dat >> Predicted.tmp; mv Predicted.tmp Predicted.dat")
write_out(structurefile, Energy_prediction(atoms, cluster_elements, support, support_size).e_total)

