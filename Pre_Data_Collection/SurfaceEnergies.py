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


'''
        CHECK for the bulk path ! ! !
'''
surf_file = "./"
try:
    surf_constrained_file = str(os.getcwd() + "/" + sys.argv[1] + "/")
except:
    surf_constrained_file = "/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Metals/Au/Surfaces/No_solvated/fcc/111/3x3/R0_L6/"
#Surf_constrained = "/home/alberto/RESEARCH/OTHER/Metals/RPBE/no_dispersion/Pure/Au/Surface/fcc/111/1x1x_Dipole/11L0R/"


def Surface_Energy(surf_constrained_file, surf_file):
    to_J = 1.60218E-19

    surf0_out = read(str(surf_constrained_file + "OUTCAR"))
    output = read(str(surf_file + "OUTCAR"))

# Checking the BULK and cell lattice agreement between slabs.
    if output.get_chemical_symbols()[0] == surf0_out.get_chemical_symbols()[0]:
        try:
            bulk_file = str("/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Metals/" + output.get_chemical_symbols()[0] + "/Bulk/fcc/")
            print("Bulk file:", bulk_file)
        except:    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            print("Non-existent ", output.get_chemical_symbols()[0], "bulk files\n")
            exit()
    else:
        print("Top surface in not proportional to bottom surface in the slab\n")
        exit()
    if sum([output.get_cell_lengths_and_angles()[i]/surf0_out.get_cell_lengths_and_angles()[i] for i in range(2)]) != 2.0:
        if sum([output.get_cell_lengths_and_angles()[i]/surf0_out.get_cell_lengths_and_angles()[i] for i in range(3, 6)]) != 3.0:
            print("Top surface in not proportional to bottom surface in the slab\n")
            print("Top surface:    ", output.get_cell_lengths_and_angles())
            print("Bottom surface: ", surf0_out.get_cell_lengths_and_angles())
            exit()

    bulk_out = read(str(bulk_file + "OUTCAR"))

    e_bulk = bulk_out.get_total_energy()			# in eV
    e_surf0 = surf0_out.get_total_energy()          # in eV
    e = output.get_total_energy()

# The formulae current employed does not requires to calculate the surface energy of the bottom (constrained) side
#    total_n_atoms, top_surf_atoms, coordination_list = Surface_Atoms(surf_constrained_file)
#    area_constrained = SurfaceArea(surf_constrained_file, top_surf_atoms)            # in m^2
#    e_surf_constrained = (e_surf0 - (len(surf0_out)/len(bulk_out)) * e_bulk) / (2*area_constrained)

    total_n_atoms, top_surf_atoms, coordination_list = Surface_Atoms(surf_file)
    area = SurfaceArea(surf_file, top_surf_atoms)				# in m^2

# area constrained and unconstrained are completely different in energy and area, e.g. pristine vs. vacancy
#    e_surf = (e - (len(output)/len(bulk_out)) * e_bulk - e_surf_constrained * area_constrained) / area
    e_surf = (e - (len(output)/len(bulk_out)) * e_bulk - (e_surf0 - (len(surf0_out)/len(bulk_out)) * e_bulk) / 2) / area

    element = output.get_chemical_symbols()[0]          # assuming that all the elements in the slab are the same

# The formulae current employed does not requires to calculate the surface energy of the bottom (constrained) side
#    return area_constrained, area, e_surf_constrained * to_J, e_surf * to_J, coordination_list, element
    return area, e_surf * to_J, coordination_list, element


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

    cutOff = neighborlist.natural_cutoffs(atoms, mult=1.25)
    a, b, d, D = neighborlist.neighbor_list('ijdD', atoms, cutOff)

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
# 4 vertex areas
            n_l = 3
            if len(j_neigh_index) < 2:
                i_2neigh_index = [b[k] for k in range(len(a)) if a[k] in i_neigh_index and b[k] != i.index]
                j_neigh_index = set([k for k in i_2neigh_index if i_2neigh_index.count(k) == 2])
                j_neigh = [k for k in range(len(a)) if a[k] == b[j] and b[k] in j_neigh_index]
                n_l = 4

            for k in sorted(j_neigh):
                done = 0                            # vertex is not done
# are i, j and k neighbours of a peak atom
                for peak in atoms_peaks:
                    peaks_neigh_index = [b[n] for n in range(len(a)) if a[n] == peak]
                    if i.index in peaks_neigh_index and b[j] in peaks_neigh_index and b[k] in peaks_neigh_index:
                        done = 1
# is this tile already measured?
                for tile in vertex_done:
                    if i.index in tile and b[j] in tile and b[k] in tile:
                        done = 1
# calculating the area and plotting the tile
                if done == 0 and n_l == 3:
                    vertex_done.append(sorted([i.index, b[j], b[k]]))
# 4 vertices tile
                if n_l == 4:
                    k_neigh = [l for l in range(len(a)) if a[l] == b[k] and b[l] in i_neigh_index and b[l] != b[j]]
                    k_neigh_index = int([b[l] for l in k_neigh][0])
# are i, j and k neighbours of a peak atom
                    for peak in atoms_peaks:
                        peaks_neigh_index = [b[n] for n in range(len(a)) if a[n] == peak]
                        if i.index in peaks_neigh_index and k_neigh_index in peaks_neigh_index and b[k] in peaks_neigh_index:
                            done = 1
                        if i.index in peaks_neigh_index and k_neigh_index in peaks_neigh_index and b[j] in peaks_neigh_index:
                            done = 1
# Bare in mind that b[k] is i.index's neighbour
                    for tile in vertex_done:
                        if i.index in tile and b[j] in tile and b[k] in tile:
                            done = 1
                        if i.index in tile and k_neigh_index in tile and b[k] in tile:
                            done = 1
                        if i.index in tile and k_neigh_index in tile and b[j] in tile:
                            done = 1
                        if b[j] in tile and k_neigh_index in tile and b[k] in tile:
                            done = 1

#                        if sorted([i.index, k_neigh_index, b[k]]) == sorted(tile):
#                            done = 1
#                        if sorted([i.index, k_neigh_index, b[j]]) == sorted(tile):
#                            done = 1
#                        if sorted([i.index, b[j], b[k]]) == sorted(tile):
#                            done = 1
#                        if sorted([k_neigh_index, b[j], b[k]]) == sorted(tile):
#                            done = 1
# calculating the area and plotting the tile
                    if done == 0 and n_l == 4:
                        vertex_done.append(sorted([i.index, k_neigh_index, b[k]]))
                        vertex_done.append(sorted([i.index, b[j], b[k]]))

    verts, area, color = Verts(atoms, vertex_done, min(z_position))
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
    while answer == "y" or type(answer) is list:
        answer = input("Would you like to remove any tile (y/n)?\n").split()
        if len(answer) < 2:                        # finished the loop when answer = n
            answer = str(answer[0])
        if answer == "y":
            vertices = input(">>> Which atoms form the tile's vertices? e.g. a b c\n").split()
            while len(vertices) != 3:
                print(">>> Only 3 vertices are accepted")
                vertices = input(">>> Which atoms form the tile's vertices?\n").split()
            else:
                vertices = sorted([int(i) for i in vertices])
                vertex_done = Remove_tile(vertices, vertex_done)
            verts, area, color = Verts(atoms, vertex_done, min(z_position))
            n_verts = len(verts)
            ax.clear
            Add_quiver_and_tiles(figure, atoms, x_max, y_max, min(z_position), a, D, color, verts)
        elif type(answer) is list and 2 < len(answer) < 4:
            vertices = sorted([int(i) for i in answer])
            vertex_done = Remove_tile(vertices, vertex_done)
            verts, area, color = Verts(atoms, vertex_done, min(z_position))
            n_verts = len(verts)
            ax.clear
            Add_quiver_and_tiles(figure, atoms, x_max, y_max, min(z_position), a, D, color, verts)

# Adding additional tiles
    answer = "y"
    while answer == "y" or type(answer) is list:
        answer = input("Would you like to cover any other area (y/n)?\n").split()
        if len(answer) < 2:                        # finished the loop when answer = n
            answer = str(answer[0])
        if answer == "y":
            vertices = input(">>> Which atoms form the tile's vertices? e.g. a b c\n").split()
            while len(vertices) != 3:
                print(">>> Only 3 vertices are accepted")
                vertices = input(">>> Which atoms form the tile's vertices?\n").split()
#            else:
            n0_verts = len(verts)
            vertex_done.append(sorted([int(i) for i in vertices]))
            verts, area, color = Verts(atoms, vertex_done, min(z_position))
            n_verts = len(verts)
            if n_verts > n0_verts:
                print("  Tile formed by", [int(i) for i in vertices], "has been added")
            else:
                vertex_done.pop()
            ax.clear
            Add_quiver_and_tiles(figure, atoms, x_max, y_max, min(z_position), a, D, color, verts)
        elif type(answer) == list and 2 < len(answer) < 4:
            n0_verts = len(verts)
            vertex_done.append(sorted([int(i) for i in answer]))
            verts, area, color = Verts(atoms, vertex_done, min(z_position))
            n_verts = len(verts)
            if n_verts > n0_verts:
                print("  Tile formed by", [int(i) for i in answer], "has been added")
            else:
                vertex_done.pop()
            ax.clear
            Add_quiver_and_tiles(figure, atoms, x_max, y_max, min(z_position), a, D, color, verts)

# Removing tiles
    answer = "y"
    while answer == "y" or type(answer) == list:
        answer = input("Would you like to remove any tile (y/n)?\n").split()
        if len(answer) < 2:                        # finished the loop when answer = n
            answer = str(answer[0])
        if answer == "y":
            vertices = input(">>> Which atoms form the tile's vertices? e.g. a b c\n").split()
            if len(vertices) != 3:
                print(">>> Only 3 vertices are accepted")
                vertices = input(">>> Which atoms form the tile's vertices?").split()
            vertices = sorted([int(i) for i in vertices])
            vertex_done = Remove_tile(vertices, vertex_done)
            verts, area, color = Verts(atoms, vertex_done, min(z_position))
            n_verts = len(verts)
            ax.clear
            Add_quiver_and_tiles(figure, atoms, x_max, y_max, min(z_position), a, D, color, verts)
        elif type(answer) == list and 2 < len(answer) < 4:
            vertices = sorted([int(i) for i in answer])
            vertex_done = Remove_tile(vertices, vertex_done)
            verts, area, color = Verts(atoms, vertex_done, min(z_position))
            n_verts = len(verts)
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

def Verts(atoms, vertex_done, z_min):
    cutOff = neighborlist.natural_cutoffs(atoms, mult=1.35)     # it has to be considerably larger for big deformations
    cutOff = sum(cutOff)/len(cutOff)*math.sqrt(2)*1.75
    z_max = max([atoms[i].position[2]-z_min for i in range(len(atoms))])
#    new_vertex = []
    verts = []
    area = []
    color = []
#    for i in range(len(vertex_done)-1):
#        new = True
#        for j in range(i+1, len(vertex_done)):
#            if sorted(vertex_done[i]) == sorted(vertex_done[j]):
#                new = False
#        if new is True:
#            new_vertex.append(sorted(vertex_done[i]))

    for tile in vertex_done:
        x = [atoms[tile[0]].position[0]]
        y = [atoms[tile[0]].position[1]]
        z = [atoms[tile[0]].position[2]-z_min]
        v = []
        d = []
        for nv in range(1, 3):
            v.append(atoms.get_distance(tile[0], tile[nv], mic=True, vector=True))
            d.append(atoms.get_distance(tile[0], tile[nv], mic=True))
            x.append(x[0]+v[nv-1][0])
            y.append(y[0]+v[nv-1][1])
            z.append(z[0]+v[nv-1][2])
#        v2 = atoms.get_distance(tile[1], tile[2], mic=True, vector=True)
#        d2 = atoms.get_distance(tile[1], tile[2], mic=True)
        if round(np.dot(v[0] / d[0], v[1] / d[1]), 5) > 1:
            print(" interatomic angle must be between pi and -pi")
        elif d[0] <= cutOff and d[1] <= cutOff:
            angle1 = np.arccos(round(np.dot(v[0] / d[0], v[1] / d[1]), 5))
#            angle2 = np.arccos(round(np.dot(v[0] / d[0], v2 / d2), 5))
#            angle3 = np.arccos(round(np.dot(v[1] / d[1], v2 / d2), 5))
#            if 0.2 < angle1 <= 2.1 and 0.2 < angle2 <= 2.1 and 0.2 < angle3 <= 2.1:     # in radiants --> 0.2 ~ 11.5 degrees; 2.1 ~ 120
#            if 0.2 < angle1 <= 2.1 and 0.2 < angle2 <= 2.1:     # in radiants --> 0.2 ~ 11.5 degrees; 2.1 ~ 120
            if 0.2 < angle1 <= 2.1:     # in radiants --> 0.2 ~ 11.5 degrees; 2.1 ~ 120
                area.append((d[0] * (d[1] * np.sin(angle1))) / 2)
                verts.append(list(zip(x, y, z)))
                if z_max > 0:
                    color.append((round(max(z), 4)-round(min(z), 4))/2/z_max)
                else:
                    color.append((round(max(z), 4)-round(min(z), 4))/2)

    return verts, area, color

def Remove_tile(vertices, vertex_done):
    new_vertex_done = []
# is this triangle already measured?
    for i, tile in enumerate(vertex_done):
        if vertices[0] in tile and vertices[1] in tile and vertices[2] in tile:
            print("  Tile formed by", vertices, "has been removed")
        else:
            new_vertex_done.append(vertex_done[i])

    return new_vertex_done

def Add_quiver_and_tiles(figure, atoms, x_max, y_max, z_min, a, D, color, verts):
    ax = figure.add_subplot(1, 1, 1, projection='3d')
# plotting the atoms, their number and the vector towards neighbours
    for i in atoms:
        ax.text(i.position[0]+x_max*0.02, i.position[1]+y_max*0.02, i.position[2]-z_min, str(i.index))
        ax.scatter3D(i.position[0], i.position[1], i.position[2]-z_min, c="k", s=50)
        i_neigh = [j for j in range(len(a)) if a[j] == i.index]
        for j in sorted(i_neigh):
            ax.quiver(i.position[0], i.position[1], i.position[2]-z_min,
                      D[j][0], D[j][1], D[j][2], length=1, arrow_length_ratio=0.05, color="b", lw=2, normalize=True)
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
# The formulae current employed does not requires to calculate the surface energy of the bottom (constrained) side
#area_constrained, area, e_surf_constrained, e_surf, coordination_list, element = Surface_Energy(surf_constrained_file, surf_file)
area, e_surf, coordination_list, element = Surface_Energy(surf_constrained_file, surf_file)

path = os.getcwd()
name = path.split("/")[-3]+"/"+path.split("/")[-2]+"/"+path.split("/")[-1]

coord_sum = []
for i in coordination_list:
    coord_sum.append(int(i)*coordination_list[i])
coord_average = sum(coord_sum)/sum([coordination_list[i] for i in coordination_list])

# The formulae current employed does not requires to calculate the surface energy of the bottom (constrained) side
#ifile = open("E_surf.dat", 'w+')
#ifile.write("# Bottom_Area(m^2)  Top_area(m^2) γ_0(J/m^2) γ_r(J/m^2) average_top_surface_coordination\n")
#ifile.write("   {:>3.9G}  {:>3.9G} {:>10.4f} {:>10.4f}   {:>3.9G}\t\t" .format(area_constrained, area, e_surf_constrained, e_surf, coord_average))
#ifile.close()

ifile = open("Trend_SEnergy.dat", 'w+')
ifile.write("# Area (Angs\u00B2)	    γ (J.m\N{SUPERSCRIPT minus}\u00B2)		Coordinations from 3 to 11	Element + path\n")
ifile.write("{:>3.9G} \t{:>10.4f} \t" .format(area*1E20, e_surf))
for i in range(3, 12):
    if str(i) in coordination_list:
        coord = coordination_list[str(i)]
    else:
        coord = 0
    ifile.write(" {:>3d}" .format(coord))
essential_surf = ["111/5x5", "111/5x5_v1", "111/5x5_v2", "111/5x5_ad1", "111/5x5_ad2", "111/5x5_ad3",
                  "001/5x5", "001/5x5_v1", "001/5x5_v2"]
if name[:7] in essential_surf:
    ifile.write("\t# {} {}\n" .format(element, name))
else:
    ifile.write("\t# {} VAL{}\n" .format(element, name))
ifile.close()
