"""
	Versions:
		Alberto: 05/2022

"""

import os
from ase.io import read
from Coordination import Coordination, Generalised_coodination
from Energies import Energies
from Properties import Areas, Cluster_surface_distance, Mean_interatomic_distance, Sphericity
from WriteData import Write_labels, write_results

""" --------------------------- CLUSTER MODEL ---------------------------"""
cluster_elements = ["Au"]  # Elements in the Cluster
isolated_cluster = "./Cluster/OUTCAR"  # file with the isolated cluster's energy
""" --------------------------- SURFACE MODEL---------------------------"""
support = "MgO"  # Support's name
support_size = [2, 2, 4]  # Dimension of the support's supercell
###################################################################################################

''' ---------------- Get and Print Results ---------------------'''
atoms = read("CONTCAR")
e_atoms = read("OUTCAR")
cluster = read(isolated_cluster)

path = os.getcwd()
name = path.split("/")[-4] + "/" + path.split("/")[-3] + "/" + path.split("/")[-2] + "/" + path.split("/")[-1]

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
          Energies(atoms, e_atoms, cluster_elements, cluster, support, support_size).e_cluster_surface,
          Energies(atoms, e_atoms, cluster_elements, cluster, support, support_size).cohesion,
          Energies(atoms, e_atoms, cluster_elements, cluster, support, support_size).adhesion,
          Energies(atoms, e_atoms, cluster_elements, cluster, support, support_size).binding,
          Energies(atoms, e_atoms, cluster_elements, cluster, support, support_size).e_total,
          name]

labels = ["N", "i_c", Coordination(atoms, cluster_elements, support).site_cluster_coordination_label, "i_cc",
          Coordination(atoms, cluster_elements, support).cluster_coord_labels, "GCN",
          Coordination(atoms, cluster_elements, support).support_cluster_min_distance_labels, "cs_dist", "cm_dist",
          "cc_dist", "shape", "sphericity", "c_i_area", "c_s_area", "Esurf", "Ecoh", "Eadh", "Eb", "Etotal",
          "  structure_path"]

Write_labels("labels.txt", labels)
write_results("data.dat", labels, values)
os.system("cat labels.txt data.dat >> Measured.dat; rm labels.txt data.dat")
