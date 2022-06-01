#!/usr/bin/env python3
"""    BASED ON Birmingham Parallel Genetic Algorithm

Please cite - A. Shayeghi et al, PCCP, 2015, 17, 2104-2112
Authors - The Johnston Group 10/7/15

"""
import os
import time
import GA.Input as In
#from GA.Graphene import Graphene                        # to generate the surface
from GA.MgO import MgO
from GA.NewPoolGA import poolGA                         # generates:  pool.dat and POSCARs
from ase.build import bulk


start_time = time.time()

""" --------------------------- MODEL ---------------------------"""
eleNames = ["Au"]                                       # elements in the cluster
eleNums = [2]                                          # number of atoms in the cluster
boxAdd = 15.0
""" SURFACE """
surfGA = True                                           # is it a supported cluster?
support = "MgO"
structure_file = "POSCAR"
surface = MgO(x=8, y=8, z=2, vac=8, clusHeight=2.0)     # how is the support's surface

""" --------------------------- ALGORITHM to generate structures ---------------------------"""
nPool = 15
cross = "random"                                        #
#cross = "weighted"                                     # determined by fitness of the two clusters selected for crossover.
# ------------------------- algorithm to generate mutants
mutType = "random"                                     # new random cluster geometry
#mutType = "move"                                        # selected from the pool and 20% of the geometry is displaced by up to 1 angstrom
#mutType = "homotop"                                    # (only bimetallics) two atoms have their atom types swapped
#mutType = "rotate"                                     # (Surface global optimisation only) selected from the pool and a random rotation is performed.
mutate = 0.4                                            # mutation ratio

r_ij = 1.2 * sum([sum(bulk(i).get_cell_lengths_and_angles()[0:3]) / 3 for i in eleNames])/len(eleNames)
eleMasses = In.masses(eleNames)
natoms = sum(eleNums)

""" --------------------------- CALCULATION ---------------------------"""
fmax = 0.5						# interatomic force maximum in the BFGS optimisation
subString = " ".join(str(i) for i in ["~/bin/RAMSCat/Predicting_MDMin.py",
            "-".join(eleNames), structure_file, support, vars(surface)['x'], vars(surface)['y'], vars(surface)['z'], fmax])      # package to calculate the Energy

StartCalc = poolGA(natoms, r_ij, eleNums, eleNames, eleMasses, mutate, nPool, cross, mutType, subString,
				   boxAdd, surface, surfGA)
# --------------------
mutType = "move"
subString = " ".join(str(i) for i in ["~/bin/RAMSCat/Predicting_MDMin.py",
            "-".join(eleNames), structure_file, support, vars(surface)['x'], vars(surface)['y'], vars(surface)['z'], fmax])      # package to calculate the Energy
StartCalc = poolGA(natoms, r_ij, eleNums, eleNames, eleMasses, mutate, nPool, cross, mutType, subString,
                                   boxAdd, surface, surfGA)
# --------------------
mutType = "rotate"
subString = " ".join(str(i) for i in ["~/bin/RAMSCat/Predicting.py",
            "-".join(eleNames), structure_file, support, vars(surface)['x'], vars(surface)['y'], vars(surface)['z'], fmax])      # package to calculate the Energy
StartCalc = poolGA(natoms, r_ij, eleNums, eleNames, eleMasses, mutate, nPool, cross, mutType, subString,
                                   boxAdd, surface, surfGA)

# --------------------

time_out = open("RAMSCat_Summary.txt", 'a+')
time_out.write("Time used to execute GA: {:>.3f} hours". format((time.time() - start_time)/3600))
time_out.write("\nThe 5 most stable structures:\n\n")
time_out.close()
os.system("grep Energy pool.dat | sort -k3n -r | head -5 >> RAMSCat_Summary.txt")
time_out = open("RAMSCat_Summary.txt", 'a+')
time_out.write("\n")
time_out.close()
os.system("head -n -1 1/Predicted.dat >> Predicted_All.dat")
os.system("re='^[0-9]+$'; for i in $(ls); do if [[ $i =~ $re ]] ; then "
		  "tail -1 $i/Predicted.dat >> Predicted_All.dat ; fi; done")
