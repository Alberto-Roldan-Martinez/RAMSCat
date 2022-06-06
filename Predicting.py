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
from Properties import Properties
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

coordination = Coordination(atoms, cluster_elements, support).coordination
generalised = Generalised_coodination(atoms, cluster_elements, coordination["Others"][0]).generalised
properties = Properties(atoms, cluster_elements, support, coordination["Others"][3], coordination["Others"][0],
                        generalised["Others"][1]).properties
energies = Energy_prediction(atoms, support, support_size, coordination["Others"][0], coordination["Others"][5],
                             coordination["Others"][1], generalised["Others"][0], generalised["Others"][1],
                             properties["c_s_area"]).energies

atoms.calc = RAMSCat(atoms, support, support_size, coordination["Others"][0], coordination["Others"][5],
                             coordination["Others"][1], generalised["Others"][0], generalised["Others"][1],
                             properties["c_s_area"])

dyn = BFGS(atoms, logfile='Optimisation.txt') #, trajectory='trajectory.traj')
#dyn = MDMin(atoms, logfile='Optimisation.txt') #, trajectory='trajectory.traj')
dyn.run(fmax=fmax, steps=500)
ase.io.vasp.write_vasp("CONTCAR.vasp", atoms, direct=False, vasp5=True, sort=True, ignore_constraints=False)

''' ---------------- Get and Print Results ---------------------'''
labels = []
values = []
for i in [coordination, generalised, properties, energies]:
    for j in list(i.keys()):
       if j != "Others":
           labels.append(j)
           values.append(i[j])

Write_labels("Predicted.tmp", labels)
write_results("Predicted.dat", labels, values)
os.system("cat Predicted.dat >> Predicted.tmp; mv Predicted.tmp Predicted.dat")
write_out(structurefile, energies["Etotal"], properties["sphericity"])
