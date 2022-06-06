"""
	Versions:
		Alberto: 05/2022

"""

import os
from ase.io import read
from Coordination import Coordination, Generalised_coodination
from Energies import Energies
from Properties import Properties
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

coordination = Coordination(atoms, cluster_elements, support).coordination
generalised = Generalised_coodination(atoms, cluster_elements, coordination["Others"][0]).generalised
properties = Properties(atoms, cluster_elements, support, coordination["Others"][3], coordination["Others"][0],
                        generalised["Others"][1]).properties
energies = Energies(atoms, e_atoms, cluster_elements, cluster, support, support_size, coordination["Others"][0],
                    generalised["Others"][1], properties["c_s_area"]).energies

labels = []
values = []
for i in [coordination, generalised, properties, energies]:
    for j in list(i.keys()):
       if j != "Others":
           labels.append(j)
           values.append(i[j])

Write_labels("labels.txt", labels)
write_results("data.dat", labels, values)
os.system("cat labels.txt data.dat >> Measured.dat; rm labels.txt data.dat")
