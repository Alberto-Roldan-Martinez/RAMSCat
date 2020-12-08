'''
    by Alberto Roldan 04/2020

'''


import os, sys, math, ase.io.vasp
import numpy as np
from ase import Atoms, neighborlist
from Library import ecoh_bulk


FileInput = "./CONTCAR"
eleNames = ["Au"]


def ClusterGCN(atoms, eleNames):
    allelements = atoms.get_chemical_symbols()
    for Cele in eleNames:
        cmask = [(Cele==ele) for ele in allelements]
    cluster = [atoms[i] for i in range(len(atoms)) if cmask[i]]               # gets the atoms in the cluster
    cluster_index = [atoms[i].index for i in range(len(atoms)) if cmask[i]]

    cutOff = neighborlist.natural_cutoffs(atoms,mult=1.2)
    a,b,d = neighborlist.neighbor_list('ijd', atoms, cutOff)

#    print(d)
    GCN = [0]*len(atoms)
    Cluster_Interface = []
    for i in range(len(a)):
        if a[i] in cluster_index:
            imask = [ (a[i]==n) for n in a]
            coordinating = [ b[i] for i in range(len(b)) if imask[i]]
            iGCN = 0

#            print ("index=",a[i],"coord=",np.bincount(a)[a[i]])

            for coord in coordinating:
                iGCN += np.bincount(a)[coord]
                if coord not in cluster_index and a[i] not in Cluster_Interface:
                    Cluster_Interface.append(a[i])
            Eb, bulk_coordination = ecoh_bulk(atoms[a[i]].symbol)
            GCN[a[i]] = iGCN/bulk_coordination
    tmp_GCN = [GCN[i] for i in range(len(GCN)) if GCN[i] != 0]
    average_GCN = sum(tmp_GCN)/len(tmp_GCN)

    return average_GCN,GCN,Cluster_Interface


atoms = ase.io.read(FileInput)
average_GCN,GCN,Cluster_Interface = ClusterGCN(atoms,eleNames)

#print(average_GCN)
