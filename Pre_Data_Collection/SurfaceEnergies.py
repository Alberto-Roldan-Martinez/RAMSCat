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
#    atoms_positions = atoms.get_scaled_positions()
#    for i in range(len(atoms)):
#        y = atoms_positions[i][1]
#        if atoms_positions[i][0] < 0:
#            atoms_positions[i][0] = atoms_positions[i][0] + 1
#        elif atoms_positions[i][0] > atoms.get_cell()[0][0] + y*np.cos(2*np.pi*atoms.get_cell_lengths_and_angles()[-1]/360):
#            atoms_positions[i][0] = atoms_positions[i][0] - 1
#        if atoms_positions[i][1] < 0:
#            atoms_positions[i][1] = atoms_positions[i][1] + 1
#        elif atoms_positions[i][1] > atoms.get_cell_lengths_and_angles()[1]:
#            atoms_positions[i][1] = atoms_positions[i][1] - 1
#
#        print(atoms[i].position, atoms_positions[i])
#    atoms.set_scaled_positions(atoms_positions)

    atoms_no_pbc = Atoms([Atom(output[i].symbol, (output[i].position)) for i in surf_atoms], cell=cell, pbc=[False])

    z_position = [output[i].position[2] for i in range(len(output)) if output[i].index in surf_atoms]
    x_max = max([output[i].position[0] for i in range(len(output)) if output[i].index in surf_atoms])
    y_max = max([output[i].position[1] for i in range(len(output)) if output[i].index in surf_atoms])

    radii = 0
    for i in surf_atoms:
        radii += sum(bulk(output[i].symbol).get_cell_lengths_and_angles()[0:3])/3
    radii = radii/len(surf_atoms)
    atoms_peaks = [atoms[i].index for i in range(len(atoms)) if atoms[i].position[2] > sum(z_position)/len(z_position)+radii*0.25]

    cutOff = neighborlist.natural_cutoffs(atoms, mult=1.2)
    a, b, d, D = neighborlist.neighbor_list('ijdD', atoms, cutOff)

    area = []
    vertex_done = []
    verts = []
    color = []
    figure = plt.figure(figsize=(10, 10), clear=True)       # prepares a figure
    ax = figure.add_subplot(1, 1, 1, projection='3d')
    for i in atoms:
        ax.text(i.position[0]+x_max*0.02, i.position[1]+y_max*0.02, i.position[2]-min(z_position), str(i.index))
        ax.scatter3D(i.position[0], i.position[1], i.position[2]-min(z_position), c="k", s=50)
        i_neigh = [j for j in range(len(a)) if a[j] == i.index]
        for j in i_neigh:
            ax.quiver(i.position[0], i.position[1], i.position[2]-min(z_position),
                      D[j][0], D[j][1], D[j][2], length=1, arrow_length_ratio=0.1, color="b", lw=2, normalize=True)

    for i in atoms:
        i_neigh = [j for j in range(len(a)) if a[j] == i.index]
        i_neigh_index = [b[j] for j in i_neigh]

        if i.index in atoms_peaks:
            print(i.index, i_neigh_index)

        for j in i_neigh:
            j_neigh_index = [b[k] for k in range(len(a)) if a[k] == b[j] and b[k] in i_neigh_index]
            j_neigh = [k for k in i_neigh if b[k] in j_neigh_index]
            for k in j_neigh:
                x = y = z = []
                if round(np.dot(D[j]/d[j], D[k]/d[k]), 5) > 1: #ij_vector/ij_distance, ik_vector/ik_distance), 5) > 1:
                    print(" interatomic angle must be between pi and -pi")
                    angle = 0
                else:
                    angle = np.arccos(round(np.dot(D[j]/d[j], D[k]/d[k]), 5)) #ij_vector/ij_distance, ik_vector/ik_distance), 5))
                if round(angle, 5) < round(np.pi, 5) and angle > 0:   # neigh -1 to have some margin
                    done = 0
                    for triangle in vertex_done:
                        if i.index in triangle and b[j] in triangle and b[k] in triangle:
                            done = 1
                    if done == 0:
                        area.append((d[j] * (d[k] * np.sin(angle))) / 2)       # N atoms contributing to Area
                        print("   triangle between  ", i.index, b[j], b[k], "  of area  ", round(area[-1], 3), "$\AA$")
                        vertex_done.append((i.index, b[j], b[k]))
                        x = [i.position[0]]
                        y = [i.position[1]]
                        z = [i.position[2]-min(z_position)]
                        x = [x[0], x[0]+D[j][0], x[0]+D[k][0]]
                        y = [y[0], y[0]+D[j][1], y[0]+D[k][1]]
                        z = [z[0], z[0]+D[j][2], z[0]+D[k][2]]
                        color.append(sum(z)/len(z))
                        verts.append(list(zip(x, y, z)))
    ax.add_collection3d(Poly3DCollection(verts, edgecolors="k", lw=0.05, facecolors=plt.cm.jet(color), alpha=0.2), zdir="z")
    ax.set_xlabel("a /$\\AA$", rotation=0, fontsize=10)
    ax.set_ylabel("b /$\\AA$", rotation=0, fontsize=10)
    ax.set_zlabel("c /$\\AA$", rotation=0, fontsize=10)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax.zaxis.set_ticklabels([])
    ax.view_init(azim=-90, elev=90)
    plt.show()

    area_total = sum(area) * 1e-20         # m^2
    print("cell area = ", output.cell[0][0] * output.cell[1][1] * 1e-20, "/$m^{2}$",
          "\ntiling area=", area_total, "/$m^{2}$", "n_area/n_verts=", len(area)/len(verts))
#    flat_area = (output.cell[0][0] * output.cell[1][1]) * 1e-20                # m^2
#    print("Area=",Area,"||","Flat_Area=",flat_area,"| Surf_atoms=",len(Surf_Atoms),"neighbouring=",neighbouring)

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
