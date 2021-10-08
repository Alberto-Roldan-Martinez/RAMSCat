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
                                                        # reads:      OUTCAR for final energy and structure
start_time = time.time()

eleNames = ["Au"]                                       # elements in the cluster
eleNums = [10]                                           # number of atoms in the cluster
boxAdd = 15.0
nPool = 15                                              # number of structures in the pool file
# --- algorithm to generate structures
cross = "random"                                        #
#cross = "weighted"										# determined by fitness of the two clusters selected for crossover.
# --- algorithm to generate mutants
#mutType = "random"										# new random cluster geometry
mutType = "move"										# selected from the pool and 20% of the geometry is displaced by up to 1 angstrom
#mutType = "homotop"									# (only bimetallics) two atoms have their atom types swapped
#mutType = "rotate"										# (Surface global optimisation only) selected from the pool and a random rotation is performed.
mutate = 0.1                                            # mutation ratio
r_ij = 2.8                                              # distance between atoms in the cluster? e.g. 3rd row -> 2.5
eleMasses = In.masses(eleNames)
natoms = sum(eleNums)

surfGA = True                                           # is it a supported cluster?
#surface = Graphene(x=8,y=8,z=1,vac=15,clusHeight=2.5)  # how is the support's surface
surface = MgO(x=8, y=8, z=4, vac=8, clusHeight=2.3)     # how is the support's surface

subString = "/home/alberto/Software/OTHER/NeuralNetwork/Predicting.py"     # package to calculate the Energy

StartCalc = poolGA(natoms, r_ij, eleNums, eleNames, eleMasses, mutate, nPool, cross, mutType, subString,
				   boxAdd, surface, surfGA)

time_out = open("RAMSCat_Summary.txt", 'w+')
time_out.write("Time used to execute GA: {:>.3f} hours". format((time.time() - start_time)/3600))
time_out.write("\nThe 5 most stable structures:\n")
time_out.close()
os.system("grep Energy pool.dat | sort -k2n | head -5 >> RAMSCat_Summary.txt")
