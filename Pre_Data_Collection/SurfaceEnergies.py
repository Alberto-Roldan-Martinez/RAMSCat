'''
    by Alberto Roldan 05/2020

'''
import os, sys, math
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
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
    surf_constrained_file = "/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Metals/Au/Surfaces/No_solvated/fcc/111/3x3/R0_L6/"
#Surf_constrained = "/home/alberto/RESEARCH/OTHER/Metals/RPBE/no_dispersion/Pure/Au/Surface/fcc/111/1x1x_Dipole/11L0R/"




def Surface_Energy(bulk_file, surf_constrained_file, surf_file):
    to_J = 1.60218E-19

    bulk_out = read(str(bulk_file + "OUTCAR"))
    surf0_out = read(str(surf_constrained_file + "OUTCAR"))
    output = read(str(surf_file + "OUTCAR"))

    e_bulk = bulk_out.get_total_energy()			# in eV
    e_surf0 = surf0_out.get_total_energy()                      # in eV
    e = output.get_total_energy()

    total_n_atoms, top_surf_atoms, coordination_list = Surface_Atoms(surf_constrained_file)
    area = SurfaceArea(surf_constrained_file, top_surf_atoms)            # in m^2
    e_surf_constrained = (e_surf0 - (len(surf0_out)/len(bulk_out)) * e_bulk) / (2*area)

    total_n_atoms, top_surf_atoms, coordination_list = Surface_Atoms(surf_file)
    area = SurfaceArea(surf_file, top_surf_atoms)				# in m^2
    e_surf = (e - (len(output)/len(bulk_out)) * e_bulk) / area - e_surf_constrained

#    print(E_bulk,E_Surf0,E)
#    if len(surf0_out) != len(output):
#        print("   Asymmetric slab: n_bulk=", len(bulk_out), "   n_S0=", len(surf0_out), "   n_S=", len(output))

    return area, e, e_surf * to_J, e_surf_constrained * to_J, (e_surf_constrained-e_surf)/e_surf_constrained*100, coordination_list


def Surface_Atoms(surf_file):
    surf_atoms_unconstrained = []
    constrained_atoms = []
    output = read(str(surf_file + "CONTCAR"))

    cutOff = neighborlist.natural_cutoffs(output, mult=1.25)
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

#    print([i for i in surf_atoms_unconstrained])

    radii = 0
    for i in surf_atoms_unconstrained:
        radii += sum(bulk(output[i].symbol).get_cell_lengths_and_angles()[0:3])/3
    radii = radii/len(surf_atoms_unconstrained)
    top_surf_atoms = [i for i in surf_atoms_unconstrained if output[i].position[2]+radii*1.2 >= sum(z_max)/len(z_max)]

#    print([i for i in top_surf_atoms])

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

#    print([i for i in surf_atoms])

    z_position = [output[i].position[2] for i in range(len(output)) if output[i].index in surf_atoms]
    x_max = max([output[i].position[0] for i in range(len(output)) if output[i].index in surf_atoms])
    y_max = max([output[i].position[1] for i in range(len(output)) if output[i].index in surf_atoms])

    cutOff = neighborlist.natural_cutoffs(atoms, mult=1.35)
    a, b, d, D = neighborlist.neighbor_list('ijdD', atoms, cutOff)

    area = []
    verts = []
    color = []
    vertex_done = []

    figure = plt.figure(figsize=(10, 10), clear=True)       # prepares a figure
    ax = figure.add_subplot(1, 1, 1, projection='3d')
# plotting the atoms, their number and the vector towards neighbours
    for i in atoms:
        ax.text(i.position[0]+x_max*0.02, i.position[1]+y_max*0.02, i.position[2]-min(z_position), str(i.index))
        ax.scatter3D(i.position[0], i.position[1], i.position[2]-min(z_position), c="k", s=50)
        i_neigh = [j for j in range(len(a)) if a[j] == i.index]
        for j in sorted(i_neigh):
            ax.quiver(i.position[0], i.position[1], i.position[2]-min(z_position),
                        D[j][0], D[j][1], D[j][2], length=1, arrow_length_ratio=0.1, color="b", lw=2, normalize=True)
# defining the peak and valley atoms : adatoms & vacancies at the surface
    radii = 0
    for i in surf_atoms:
        radii += sum(bulk(output[i].symbol).get_cell_lengths_and_angles()[0:3])/3
    radii = radii/len(surf_atoms)
    atoms_peaks = [atoms[i].index for i in range(len(atoms)) if atoms[i].position[2] > sum(z_position)/len(z_position)+radii*0.25]
    atoms_valleys = [atoms[i].index for i in range(len(atoms)) if atoms[i].position[2] < sum(z_position)/len(z_position)-radii*0.25]
# finding neighbours of i and neighbours of neighbours that are also neighbours of i
    for i in atoms:
        i_neigh = [j for j in range(len(a)) if a[j] == i.index]
        i_neigh_index = [b[j] for j in i_neigh]
        for j in sorted(i_neigh):
            j_neigh_index = [b[k] for k in range(len(a)) if a[k] == b[j] and b[k] in i_neigh_index]
            j_neigh = [k for k in i_neigh if b[k] in j_neigh_index]
# are i, j and k neighbours of a valley atom
            if i.index in atoms_valleys and b[j] in atoms_valleys:
                j_exclusive = []
                j_exclusive_neigh = [k for k in range(len(a)) if a[k] == b[j] and b[k] not in i_neigh_index and b[k] not in atoms_valleys]
                for k in j_exclusive_neigh:
                    k_neigh = [n for n in range(len(a)) if a[n] == b[k] if b[n] in i_neigh_index and b[n] not in atoms_valleys]
                    if len(k_neigh) > 0:
                        j_exclusive.append(k)
                ik_distance = []
                for k in j_exclusive:
                    ik_distance.append(atoms.get_distance(i.index, b[k], mic=True))
                for k in j_exclusive:
                    dist = atoms.get_distance(i.index, b[k], mic=True)
                    if dist == min(ik_distance) and dist < sum(cutOff)/len(cutOff)*2.75:
                        j_neigh_index.append(b[k])
                        j_neigh.append(k)
            for k in sorted(j_neigh):
                x = y = z = []
                done = 0                            # vertex is not done
# are i, j and k neighbours of a peak atom
                for peak in atoms_peaks:
                    peaks_neigh_index = [b[n] for n in range(len(a)) if a[n] == peak]
                    if i.index in peaks_neigh_index and b[j] in peaks_neigh_index and b[k] in peaks_neigh_index:
                        done = 1
# is this triangle already measured?
                for triangle in vertex_done:
                    if i.index in triangle and b[j] in triangle and b[k] in triangle:
                        done = 1
# calculating the area and plotting the triangle
                if done == 0:
                    area, color, verts, vertex_done = Add_triangle(atoms, [i.index, b[j], b[k]], min(z_position), max(z_position),
                                                      color, verts, area, vertex_done)
    n_verts = len(verts)
    ax.add_collection3d(Poly3DCollection(verts, edgecolors="k", lw=0.1, facecolors=plt.cm.inferno(color), alpha=0.4), zdir="z")
#    ax.set_xlabel("a /$\\AA$", rotation=0, fontsize=10); ax.set_ylabel("b /$\\AA$", rotation=0, fontsize=10); ax.set_zlabel("c /$\\AA$", rotation=0, fontsize=10)
#    ax.xaxis.set_ticklabels([]); ax.yaxis.set_ticklabels([]); ax.zaxis.set_ticklabels([])
    figure.patch.set_visible(False)
    ax.axis('off')
    ax.view_init(azim=-90, elev=90)
    plt.ion()
    ax.set_xlim3d([0, x_max])
    ax.set_ylim3d([0, y_max])
    ax.set_zlim3d([0, (x_max+y_max)/2])
    plt.show()

# Removing tiles
    answer = "y"
    while answer == "y":
        answer = str(input("Would you like to remove any tile (y/n)?\n"))
        if answer == "y":
            vertices = input("   Which three atoms form the tile's vertices? e.g. a b c\n").split()
            vertices = [int(i) for i in vertices]
            if len(vertices) != 3:
                print(">>> Only three vertices are accepted")
                vertices = input(">>> Which three atoms form the tile's vertices?").split()
                vertices = [int(i) for i in vertices]
            area, color, verts, vertex_done = Remove_tile(vertices, vertex_done, color, verts, area)
            n_verts -= 1
            ax.clear
            Add_quiver_and_tiles(figure, atoms, x_max, y_max, min(z_position), a, D, color, verts)

# Adding additional tiles
    answer = "y"
    while answer == "y":
        answer = str(input("Would you like to cover any other area (y/n)?\n"))
        if answer == "y":
            vertices = input("   Which three atoms form the tile's vertices? e.g. a b c\n").split()
            vertices = [int(i) for i in vertices]
            if len(vertices) != 3:
                print(">>> Only three vertices are accepted")
                vertices = input(">>> Which three atoms form the tile's vertices?").split()
                vertices = [int(i) for i in vertices]
            area, color, verts, vertex_done = Add_triangle(atoms, vertices, min(z_position), max(z_position), color, verts, area, vertex_done)
            n_verts += 1
            ax.clear
            Add_quiver_and_tiles(figure, atoms, x_max, y_max, min(z_position), a, D, color, verts)

# Removing tiles
    answer = "y"
    while answer == "y":
        answer = str(input("Would you like to remove any tile (y/n)?\n"))
        if answer == "y":
            vertices = input("   Which three atoms form the tile's vertices? e.g. a b c\n").split()
            vertices = [int(i) for i in vertices]
            if len(vertices) != 3:
                print(">>> Only three vertices are accepted")
                vertices = input(">>> Which three atoms form the tile's vertices?").split()
                vertices = [int(i) for i in vertices]
            area, color, verts, vertex_done = Remove_tile(vertices, vertex_done, color, verts, area)
            n_verts -= 1
            ax.clear
            Add_quiver_and_tiles(figure, atoms, x_max, y_max, min(z_position), a, D, color, verts)

    if len(area)/n_verts != 1:
        print("n_area/n_verts=", len(area)/n_verts, "\n   n areas=", len(area), "\n   n vertex=", n_verts)
    area_total = sum(area) * 1e-20         # m^2
    print("cell area = ", output.cell[0][0] * output.cell[1][1] * 1e-20, "/$m^{2}$",
          "\ntiling area=", area_total, "/$m^{2}$\n--------------------------------------------\n")
#    flat_area = (output.cell[0][0] * output.cell[1][1]) * 1e-20                # m^2
#    print("Area=",Area,"||","Flat_Area=",flat_area,"| Surf_atoms=",len(Surf_Atoms),"neighbouring=",neighbouring)

    return area_total

def Add_triangle(atoms, vertices, z_min, z_max, color, verts, area, vertex_done):
    x = [atoms[vertices[0]].position[0]]
    y = [atoms[vertices[0]].position[1]]
    z = [atoms[vertices[0]].position[2]-z_min]
    v1 = atoms.get_distance(vertices[0], vertices[1], mic=True, vector=True)
    d1 = atoms.get_distance(vertices[0], vertices[1], mic=True)
    v2 = atoms.get_distance(vertices[0], vertices[2], mic=True, vector=True)
    d2 = atoms.get_distance(vertices[0], vertices[2], mic=True)
    x = [x[0], x[0]+v1[0], x[0]+v2[0]]
    y = [y[0], y[0]+v1[1], y[0]+v2[1]]
    z = [z[0], z[0]+v1[2], z[0]+v2[2]]
    if round(np.dot(v1 / d1, v2 / d2), 5) > 1:
        print(" interatomic angle must be between pi and -pi")
        angle = 0
    else:
        angle = np.arccos(round(np.dot(v1 / d1, v2 / d2), 5))
    if round(angle, 5) <= round(np.pi/1.25, 5) and angle > 1E-3:
        area.append((d1 * (d2 * np.sin(angle))) / 2)       # N atoms contributing to Area
        verts.append(list(zip(x, y, z)))
        vertex_done.append(vertices)
        color.append((sum(z)/len(z))/((z_max-z_min)*2))

    return area, color, verts, vertex_done

def Remove_tile(vertices, vertex_done, color, verts, area):
    new_verts = []
    new_color = []
    new_area = []
    new_vertex_done = []
# is this triangle already measured?
    for i, tile in enumerate(vertex_done):
        if vertices[0] in tile and vertices[1] in tile and vertices[2] in tile:
            print("  Tile formed by", vertices, "has been removed")
        else:
            new_verts.append(verts[i])
            new_color.append(color[i])
            new_area.append(area[i])
            new_vertex_done.append(vertex_done[i])

    return new_area, new_color, new_verts, new_vertex_done

def Add_quiver_and_tiles(figure, atoms, x_max, y_max, z_min, a, D, color, verts):
    ax = figure.add_subplot(1, 1, 1, projection='3d')
# plotting the atoms, their number and the vector towards neighbours
    for i in atoms:
        ax.text(i.position[0]+x_max*0.02, i.position[1]+y_max*0.02, i.position[2]-z_min, str(i.index))
        ax.scatter3D(i.position[0], i.position[1], i.position[2]-z_min, c="k", s=50)
        i_neigh = [j for j in range(len(a)) if a[j] == i.index]
        for j in sorted(i_neigh):
            ax.quiver(i.position[0], i.position[1], i.position[2]-z_min,
                      D[j][0], D[j][1], D[j][2], length=1, arrow_length_ratio=0.1, color="b", lw=2, normalize=True)
    ax.add_collection3d(Poly3DCollection(verts, edgecolors="k", lw=0.1, facecolors=plt.cm.inferno(color), alpha=0.4), zdir="z")
    figure.patch.set_visible(False)
    ax.axis('off')
    ax.view_init(azim=-90, elev=90)
#    plt.ion()
    ax.set_xlim3d([0, x_max])
    ax.set_ylim3d([0, y_max])
    ax.set_zlim3d([0, (x_max+y_max)/2])
    plt.show()



##############################################################################################

area, e, e_surf, e_surf_constrained, relaxation, coordination_list = Surface_Energy(bulk_file, surf_constrained_file, surf_file)
ifile = open("E_surf.dat", 'w+')
#ifile.write("# Tatoms Satoms Area(m^2) E(eV/atom) γ_0(J/m^2) γ_r(J/m^2) (γ_0-γ_r)/γ_0(%)  || Coordination_list\n")
#ifile.write("    %d      %d     %.3G    %.3f      %.3f      %.3f         %.3f   " % (total_n_atoms, len(top_surf_atoms),
#                                                                                     area, e/total_n_atoms,
#                                                                                     e_surf_constrained, e_surf,
#                                                                                     relaxation))
ifile.write("# Area(m^2)   γ_r(J/m^2) || Coordination_list\n")
ifile.write("  %.5G       %.3f " % (area, e_surf))
ifile.write("||  ")
for i in coordination_list:
    ifile.write("%5s x %d" %(i, coordination_list[i]))
ifile.write("\n")
ifile.close()
