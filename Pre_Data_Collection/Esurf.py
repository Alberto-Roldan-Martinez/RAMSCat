'''
    by Alberto Roldan 05/2020

'''
import os, sys, math
import numpy as np
from array import *
from ase.io import read, write
from ase import Atoms, Atom, neighborlist
from ase.build import bulk
from Coordination import GCN


Surf = "./OUTCAR"
Bulk = "/home/alberto/RESEARCH/OTHER/Metals/rPBE/no_dispersion/Pure/Au/Bulk/fcc/K01/OUTCAR"
try:
    Surf_constrained = str(os.getcwd() + "/" + sys.argv[1] + "/OUTCAR")
except:
    Surf_constrained = "/home/alberto/RESEARCH/OTHER/Metals/rPBE/no_dispersion/Pure/Au/Surface/fcc/111/1x1/Layers/Vacuum/K02/0R/15L/OUTCAR"
#Surf_constrained = "/home/alberto/RESEARCH/OTHER/Metals/rPBE/no_dispersion/Pure/Au/Surface/fcc/111/1x1x_Dipole/11L0R/OUTCAR"




def Surface_Energy(Bulk, Surf_constrained, Surf, Surf_Atoms):
    to_J = 1.60218E-19

    bulk_out = read(Bulk)
    Surf0_out = read(Surf_constrained)
    output = read(Surf)

    E_bulk = bulk_out.get_total_energy()			# in eV
    E_Surf0 = Surf0_out.get_total_energy()                      # in eV
    E = output.get_total_energy()

    Area = SurfaceArea(Surf_constrained, Surf_Atoms)            # in m^2
    E_surf_constrained = (E_Surf0 - (len(Surf0_out)/len(bulk_out)) * E_bulk) / (2*Area)

    Area = SurfaceArea(Surf, Surf_Atoms)				# in m^2
    E_surf = (E - (len(output)/len(bulk_out)) * E_bulk) / Area - E_surf_constrained

#    print(E_bulk,E_Surf0,E)
    if len(Surf0_out) != len(output):
        print("   n_bulk=", len(bulk_out), "   n_S0=", len(Surf0_out), "   n_S=", len(output))

    return Area, E, E_surf * to_J, E_surf_constrained * to_J, (E_surf_constrained-E_surf)/E_surf_constrained*100


def Surface_Atoms(Surf):
    Surf_Atoms_unconstrained = []
    constrained_atoms = []
    output = read(Surf)
    iGCN = GCN(Surf, output.get_chemical_symbols())
    coord = iGCN.coordination_list
    gcn_sum = 0
    gcn_list = {}
    coordination_list = {}

    constrainment = [i.index for i in output.constraints]
    if len(constrainment) != 0:
        for i in constrainment[0]:
            constrained_atoms.append(i)

    Zsum = []
    for atom in output:
        if coord[str(atom.index)] < 12 and atom.index not in constrained_atoms:
            Surf_Atoms_unconstrained.append(atom.index)
            Zsum.append(atom.position[2])

    Surf_Atoms = [i for i in Surf_Atoms_unconstrained if output[i].position[2]+2 >= sum(Zsum)/len(Zsum)]
#    Surf_Atoms = Surf_Atoms_unconstrained

    for i in Surf_Atoms:
        gcn_sum += iGCN.gcn[atom.index]

    if len(Surf_Atoms) == 0:
        Surf_Atoms.append(output[-1].index)
        gcn_sum += iGCN.gcn[output[-1].index]

    surf_GCN = gcn_sum/len(Surf_Atoms)

    for index in Surf_Atoms:
        if str(round(iGCN.gcn[index],3)) not in gcn_list:
            gcn_list[str(round(iGCN.gcn[index], 3))] = 1
        else:
            gcn_list[str(round(iGCN.gcn[index], 3))] += 1

        if str(coord[str(index)]) not in coordination_list:
            coordination_list[str(coord[str(index)])] = 1
        else:
            coordination_list[str(coord[str(index)])] += 1

    sorted(gcn_list)
    sorted(coordination_list)

    return len(output), Surf_Atoms, surf_GCN, gcn_list, coordination_list


def SurfaceArea(Surf, Surf_Atoms):

    output = read(Surf)
    atoms = Atoms([Atom(output[i].symbol, (output[i].position)) for i in Surf_Atoms])
    atoms.set_cell(output.get_cell())
    atoms.set_pbc(output.get_pbc())

    Z_position = [i.position[2] for i in atoms]
    radii = sum(sum([bulk(i.symbol).get_cell_lengths_and_angles()[0:3] for i in atoms])/3)/len(atoms)

    if max(Z_position) - min(Z_position) < radii*0.6:
        print("Smooth surface! (", round(radii, 3), "*60% <", round(max(Z_position) - min(Z_position), 3), "Angs)")
        Area = (output.cell[0][0] * output.cell[1][1]) * 1e-20               # m^2
    else:
        print("Rough surface! (", round(radii, 3), "*60% >", round(max(Z_position) - min(Z_position), 3), "Angs)")
        cutOff = neighborlist.natural_cutoffs(atoms, mult=1.2)
        a, b, d, D =  neighborlist.neighbor_list('ijdD', atoms, cutOff)

        triangular_area = {}
        Area = []
        neighbouring = 0 
        for i in atoms:
            neigh = [j for j in range(len(a)) if a[j] == i.index]
            neighbouring += len(neigh)/len(atoms)
            for j in neigh:
                for k in neigh:
                    if round(np.dot(D[j]/d[j], D[k]/d[k]), 5) > 1:
                        print(" interatomic angle must be between pi and -pi")
                        angle = 0
                    else:
                        angle = np.arccos(round(np.dot(D[j]/d[j], D[k]/d[k]), 5))
                    if round(angle, 5) <= round(2*np.pi/(len(neigh)-1), 5) and angle > 1E-3:   # neigh -1 to have some margin
                        heigh = d[k] * np.sin(angle)
                        Area.append(((d[j] * heigh) / 2)/len(neigh))       # N atoms contributing to Area
                        key = "i" + str(i.index) + "j" + str(j) + "k" + str(a[k])
                        triangular_area.setdefault(key, []).append(Area[-1])

        Area = sum(Area) * 1e-20         # m^2
#        flat_area = (output.cell[0][0] * output.cell[1][1]) * 1e-20                # m^2
#        print("Area=",Area,"||","Flat_Area=",flat_area,"| Surf_atoms=",len(Surf_Atoms),"neighbouring=",neighbouring)


    return Area 


##############################################################################################


Tatoms, Surf_Atoms, surf_GCN, gcn_list, coordination_list = Surface_Atoms(Surf)
Area, E, E_surf, E_surf_constrained, Relaxation = Surface_Energy(Bulk, Surf_constrained, Surf, Surf_Atoms)


ifile = open("E_surf.dat", 'w+')
ifile.write("# Tatoms Satoms ave_GCN Area(m^2) E(eV/atom) γ_0(J/m^2) γ_r(J/m^2) (γ_0-γ_r)/γ_0(%)      gcn_list  || Coordination_list\n")
ifile.write("    %d      %d     %.3f   %.3G    %.3f      %.3f      %.3f         %.3f   " 	%(Tatoms, len(Surf_Atoms), surf_GCN, Area, E/Tatoms, E_surf_constrained, E_surf, Relaxation))
for i in gcn_list:
    ifile.write("  %7s x %d" %(i,gcn_list[i]))
ifile.write("   ||")
for i in coordination_list:
    ifile.write("  %7s x %d" %(i, coordination_list[i]))
ifile.write("\n")
ifile.close()
