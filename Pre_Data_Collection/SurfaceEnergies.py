'''
    by Alberto Roldan 05/2020

'''
import os, sys, math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from ase.io import read #, write
from ase import Atoms, Atom, neighborlist
from ase.build import bulk



surf_file = "./"
bulk_file = "/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Metals/Au/Bulk/fcc/"
try:
    surf_constrained_file = str(os.getcwd() + "/" + sys.argv[1] + "/")
except:
    surf_constrained_file = "/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Metals/Au/Surfaces/No_solvated/fcc/111/2x2/R0_L6/"
#Surf_constrained = "/home/alberto/RESEARCH/OTHER/Metals/RPBE/no_dispersion/Pure/Au/Surface/fcc/111/1x1x_Dipole/11L0R/"




def Surface_Energy(bulk_file, surf_constrained_file, surf_file, top_surf_atoms):
    to_J = 1.60218E-19

    bulk_out = read(str(bulk_file + "OUTCAR"))
    surf0_out = read(str(surf_constrained_file + "OUTCAR"))
    output = read(str(surf_file + "OUTCAR"))

    e_bulk = bulk_out.get_total_energy()			# in eV
    e_surf0 = surf0_out.get_total_energy()                      # in eV
    e = output.get_total_energy()

    area = SurfaceArea(surf_constrained_file, top_surf_atoms)            # in m^2
    e_surf_constrained = (e_surf0 - (len(surf0_out)/len(bulk_out)) * e_bulk) / (2*area)

    area = SurfaceArea(surf_file, top_surf_atoms)				# in m^2
    e_surf = (e - (len(output)/len(bulk_out)) * e_bulk) / area - e_surf_constrained

#    print(E_bulk,E_Surf0,E)
    if len(surf0_out) != len(output):
        print("   n_bulk=", len(bulk_out), "   n_S0=", len(surf0_out), "   n_S=", len(output))

    return area, e, e_surf * to_J, e_surf_constrained * to_J, (e_surf_constrained-e_surf)/e_surf_constrained*100


def Surface_Atoms(surf_file):
    surf_atoms_unconstrained = []
    surf_atoms_constrained = []
    constrained_atoms = []
    output = read(str(surf_file + "CONTCAR"))

    cutOff = neighborlist.natural_cutoffs(output, mult=1.2)
    a,b,d = neighborlist.neighbor_list('ijd', output, cutOff)


    constrainment = [i.index for i in output.constraints]
    if len(constrainment) != 0:
        for i in constrainment[0]:
            constrained_atoms.append(i)

    z_max = []
    for atom in output:
        if np.bincount(a)[atom.index] < 12 and atom.index not in constrainment[0]:
                surf_atoms_unconstrained.append(atom.index)
                z_max.append(atom.position[2])

    top_surf_atoms = [i for i in surf_atoms_unconstrained if output[i].position[2]+2 >= sum(z_max)/len(z_max)]

    coordination_list = {}
    for index in top_surf_atoms:
        if str(np.bincount(a)[index]) not in coordination_list:
            coordination_list[str(np.bincount(a)[index])] = 1
        else:
            coordination_list[str(np.bincount(a)[index])] += 1
    sorted(coordination_list)

    return len(output), top_surf_atoms, coordination_list


def SurfaceArea(slab_file, surf_atoms):
    output = read(str(slab_file + "OUTCAR"))
    cell = output.get_cell()
    cell_pbc = output.get_pbc()
    output = read(str(slab_file + "CONTCAR"))
    atoms = Atoms([Atom(output[i].symbol, (output[i].position)) for i in surf_atoms], cell=cell, pbc=cell_pbc)
    atoms_no_pbc = Atoms([Atom(output[i].symbol, (output[i].position)) for i in surf_atoms], cell=cell, pbc=[False])

    z_position = [output[i].position[2] for i in range(len(output)) if output[i].index in surf_atoms]
    z_max = max(z_position)-min(z_position)
    x_max = max([output[i].position[0] for i in range(len(output)) if output[i].index in surf_atoms])
    y_max = max([output[i].position[1] for i in range(len(output)) if output[i].index in surf_atoms])

    radii = 0
    for i in surf_atoms:
        radii += sum(bulk(output[i].symbol).get_cell_lengths_and_angles()[0:3])/3
    radii = radii/len(surf_atoms)

    if max(z_position) - min(z_position) < radii*0.5:
        print("Smooth surface! (", round(radii, 3), "*50% <", round(max(z_position) - min(z_position), 3), "Angs)")
        area_total = (output.cell[0][0] * output.cell[1][1]) * 1e-20               # m^2
    else:
        print("Rough surface! (", round(radii, 3), "*50% >=", round(max(z_position) - min(z_position), 3), "Angs)")
        cutOff = neighborlist.natural_cutoffs(atoms, mult=1.6)
        a, b = neighborlist.neighbor_list('ij', atoms, cutOff)
        a_no_pbc, b_no_pbc = neighborlist.neighbor_list('ij', atoms_no_pbc, cutOff)

        area = []
        vertex_done = []
        verts = []
        verts_extra = []
        color = []
        figure = plt.figure(figsize=(10, 10), clear=True)       # prepares a figure
        ax = figure.add_subplot(1, 1, 1, projection='3d')
        for i in atoms:
            i_neigh = [b[j] for j in range(len(a)) if a[j] == i.index]
            neighbour_no_pbc = [b_no_pbc[j] for j in range(len(a_no_pbc)) if a_no_pbc[j] == i.index]
            x = []
            y = []
            z = []
            for j in i_neigh:
                j_neigh = [b[k] for k in range(len(a)) if a[k] == atoms[j].index and b[k] in i_neigh]
                for k in j_neigh:
                    ij_distance = atoms.get_distance(i.index, j)
                    ij_vector = atoms.get_distance(i.index, j, vector=True)
                    ik_distance = atoms.get_distance(i.index, k)
                    ik_vector = atoms.get_distance(i.index, k, vector=True)
                    if round(np.dot(ij_vector/ij_distance, ik_vector/ik_distance), 5) > 1:
                        print(" interatomic angle must be between pi and -pi")
                        angle = 0
                    else:
                        angle = np.arccos(round(np.dot(ij_vector/ij_distance, ik_vector/ik_distance), 5))
                    if round(angle, 5) <= round(2*np.pi/(len(i_neigh)-1), 5) and angle > 1E-3:   # neigh -1 to have some margin
                        ii = i.index
                        jj = atoms[j].index
                        kk = atoms[k].index
                        if str(ii)+str(jj)+str(kk) not in vertex_done:
                            area.append(((ij_distance * (ik_distance * np.sin(angle))) / 2)/len(i_neigh))       # N atoms contributing to Area

                            vertex_done.append(str(ii)+str(jj)+str(kk))
                            vertex_done.append(str(ii)+str(kk)+str(jj))
                            vertex_done.append(str(jj)+str(kk)+str(ii))
                            vertex_done.append(str(jj)+str(ii)+str(kk))
                            vertex_done.append(str(kk)+str(jj)+str(ii))
                            vertex_done.append(str(kk)+str(ii)+str(jj))
                            if j in neighbour_no_pbc and k in neighbour_no_pbc:
                                x = [i.position[0], atoms[j].position[0], atoms[k].position[0]]
                                y = [i.position[1], atoms[j].position[1], atoms[k].position[1]]
                                z = [i.position[2], atoms[j].position[2], atoms[k].position[2]]
                                x = [round(value/x_max, 2) for value in x]
                                y = [round(value/y_max, 2) for value in y]
                                z = [round((value-min(z_position))/z_max, 2) for value in z]
                                color.append(sum(z)/len(z))
                                verts.append(list(zip(x, y, z)))
                                ax.add_collection3d(Poly3DCollection(verts, edgecolors="k", facecolors=plt.cm.jet(color)), zdir="z")

                            x2 = [round(value/x_max + 1, 2) for value in x]
#                            verts_extra.append(list(zip(x2, y, z)))
                            y2 = [round(value/y_max + 1, 2) for value in y]
#                            verts_extra.append(list(zip(x, y2, z)))
#                            verts_extra.append(list(zip(x2, y2, z)))
#                            ax.add_collection3d(Poly3DCollection(verts_extra, edgecolors="k", facecolors=plt.cm.jet(color), alpha=0.5), zdir="z")
#                        else:
#                            print(str(ii)+str(jj)+str(kk)," is done!")
        print("area=",len(area),len(verts))
        x = [i.position[0]/x_max for i in atoms]
        y = [i.position[1]/y_max for i in atoms]
        z = [(i.position[2]-min(z_position))/z_max for i in atoms]
        ax.scatter3D(x, y, z, c="k", s=50)#, edgecolor='none', linewidth=0, antialiased=False)
        ax.set_xlabel("a /$\\AA$", rotation=0, fontsize=10)
        ax.set_ylabel("b /$\\AA$", rotation=0, fontsize=10)
        ax.set_zlabel("c /$\\AA$", rotation=0, fontsize=10)
#        ax.xaxis.set_ticklabels([])
#        ax.yaxis.set_ticklabels([])
        ax.zaxis.set_ticklabels([])
        ax.view_init(azim=-90, elev=90)
        plt.show()

        area_total = sum(area) * 1e-20         # m^2
#        flat_area = (output.cell[0][0] * output.cell[1][1]) * 1e-20                # m^2
#        print("Area=",Area,"||","Flat_Area=",flat_area,"| Surf_atoms=",len(Surf_Atoms),"neighbouring=",neighbouring)

    return area_total


##############################################################################################

total_n_atoms, top_surf_atoms, coordination_list = Surface_Atoms(surf_file)
area, e, e_surf, e_surf_constrained, relaxation = Surface_Energy(bulk_file, surf_constrained_file, surf_file,
                                                                 top_surf_atoms)
ifile = open("E_surf.dat", 'w+')
#ifile.write("# Tatoms Satoms Area(m^2) E(eV/atom) γ_0(J/m^2) γ_r(J/m^2) (γ_0-γ_r)/γ_0(%)  || Coordination_list\n")
#ifile.write("    %d      %d     %.3G    %.3f      %.3f      %.3f         %.3f   " % (total_n_atoms, len(top_surf_atoms),
#                                                                                     area, e/total_n_atoms,
#                                                                                     e_surf_constrained, e_surf,
#                                                                                     relaxation))
ifile.write("# Area(m^2) γ_r(J/m^2) || Coordination_list\n")
ifile.write("  %.3G       %.3f " % (area, e_surf))
ifile.write("||  ")
for i in coordination_list:
    ifile.write("%5s x %d" %(i, coordination_list[i]))
ifile.write("\n")
ifile.close()
