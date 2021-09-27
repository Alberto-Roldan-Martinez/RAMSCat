'''    BASED ON Birmingham Parallel Genetic Algorithm

Please cite - A. Shayeghi et al, PCCP, 2015, 17, 2104-2112
Authors - The Johnston Group 10/7/15

'''
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
nPool = 500                                              # number of structures in the pool file
cross = "random"                                        # algorithm to generate structures
mutType = "move"                                        # algorithm to generate mutants
mutate = 0.1                                            # mutation ratio
r_ij = 2.8                                              # distance between atoms in the cluster?
eleMasses = In.masses(eleNames)
natoms = sum(eleNums)

surfGA = True                                           # is it a supported cluster?
#surface = Graphene(x=8,y=8,z=1,vac=15,clusHeight=2.5)  # how is the support's surface
surface = MgO(x=8, y=8, z=2, vac=6, clusHeight=2.3)     # how is the support's surface

subString = "/home/alberto/Software/OTHER/NeuralNetwork/Predicting.py"     # package to calculate the Energy

StartCalc = poolGA(natoms, r_ij, eleNums, eleNames, eleMasses, mutate, nPool, cross, mutType, subString,
				   boxAdd, surface, surfGA)

time_out = open("Calculation_time.txt", 'w+')
time_out.write("Time used to execute GA: {} seconds". format(time.time() - start_time))
time_out.close()