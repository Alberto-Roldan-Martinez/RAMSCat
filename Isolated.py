

''' the 3 is coord += math.exp(-(3*Ndistance)**2) 
    is discretional. It should be exchanged by a 
    function describing the elasticity  '''



import os, sys, math, ase.io.vasp
import numpy as np
from ase import Atoms, neighborlist
from ase.io import read
from ase.build import bulk,molecule
#from ase.data import atomic_numbers, covalent_radii
from Library import OptAtomdistance,IsolatedAtoms
from WriteData import WriteLabels,WriteResults

isolatedSystem = "OUTCAR"            # ".../OTHER/Supported/MgO/Au/Basin_Hopping/1Au/gas/OUTCAR"
isolatedCluster = "CONTCAR"

system = ase.io.vasp.read_vasp_out(isolatedSystem, index=-1)
atoms = ase.io.read(isolatedCluster)

Ecluster = system.get_total_energy()
allelements = system.get_chemical_symbols()

cluster = atoms

coord=0
Eatoms=0
eleNames=[allelements[0]]
new = 1
for ele in allelements:
    for nam in eleNames:
        if ele == nam:
            new = 0
    if new == 1:
        eleNames.append(ele)

for Cele in eleNames:
    Eatoms += IsolatedAtoms(Cele)

cutOff = neighborlist.natural_cutoffs(atoms, mult=1.2) # = [ (covalent_radii[atom.number] * 1.2) for atom in atoms]
coord = neighborlist.neighbor_list('i', atoms, cutOff)

cc = sum(np.bincount(coord))/len(np.bincount(coord))
clusterSize = len(cluster)
Ecohesion = (Ecluster - Eatoms)/len(cluster)

columnLabels=["N","cc", "Ecoh(eV/atom)", "Etotal(eV)"]
values=[clusterSize,cc,Ecohesion,Ecluster]

WriteLabels("gaslabels.txt",columnLabels)
WriteResults("gasdata.dat",values)
