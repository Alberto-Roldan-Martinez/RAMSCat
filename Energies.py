"""
    by Alberto Roldan

    Reads inputs containing the energies and calculate
        Ecohesion
        Eadhesion
        TotalEnergy

    Versions:
        Alberto: 08/2019


    STRUCTURE:

"""

import numpy as np
from ase.io import read
from Coordination import Coordination
from GCN import Generalised_coodination
from Areas import Areas
from Library import isolated_atoms, surf_energies, ecoh_trend, e_adh_energies, supports


class Energies:
    def __init__(self, inputfiles, isolated_support, isolated_cluster, cluster_elements, support):
        system = read(inputfiles[0], index=-1)
        cluster = read(isolated_cluster)
        slab = read(isolated_support)

        e_system = system.get_total_energy()
        e_cluster = cluster.get_total_energy()
        e_slab = slab.get_total_energy()
        e_atoms = sum([isolated_atoms(system[i].symbol) for i in range(len(system))
                       if system[i].symbol in cluster_elements])

# c_coord = dictionary with the indexes of coordinating atoms within the cluster
        coordination = Coordination(inputfiles[1], cluster_elements, support)
        c_coord = coordination.cluster_coordinating
# c_surf = indexes of cluster atoms with coordination within the cluster lower than its bulk
        gcn = Generalised_coodination(inputfiles[1], cluster_elements, support)
        c_surf = gcn.cluster_surface_index
# c_surf_area = area exposed to the vacuum according to interpolated area/atom in the Library
        area = Areas(inputfiles[1], cluster_elements, support)
        c_surf_area = area.cluster_surface_area

# cohesion = cohesion energy per atom in the cluster in eV/atom from the DFT data
# adhesion = adhesion energy in eV from the DFT calculated data
# binding = binding energy in eV from the DFT calculated data
# e_total = total energy in eV of the supported cluster from the DFT data
        self.cohesion = float((e_cluster - e_atoms) / len(cluster))
        self.adhesion = float(e_system - (e_cluster + e_slab))
#        if sum(coordination.site_cluster_coordination) > 0:
#            self.normalised_adhesion = self.adhesion / sum(coordination.site_cluster_coordination)
#        else:
#            self.normalised_adhesion = self.adhesion
        self.binding = float((e_system - (e_slab + e_atoms))/len(cluster))
        self.e_total = float(e_system)
# e_cluster_surface = surface energy (SE) in J/m^2 of the cluster atoms expossed to vacuum according to interpolated SE and area / atom in the Library
        self.e_cluster_surface = surface_energy(system, c_coord, c_surf, c_surf_area)


def surface_energy(system, c_coord, c_surf, c_surf_area):
    ev_to_joules = 1.60218E-19
    e_c_surf = sum([surf_energies(system[i].symbol, len(c_coord[str(i)])) for i in c_surf])
    return float(e_c_surf * ev_to_joules / (c_surf_area * 1e-20))                  # J/m^2


class Energy_prediction:
    def __init__(self, inputfile, cluster_elements, support, support_size):
        system = read(inputfile)
        e_slab = supports(support, support_size)

# c_coord = dictionary with the indexes of coordinating atoms within the cluster
# interface_distances = dictionary of interface cluster atoms with the minimum distances to site X and Y
# interface_indexes = dictionary with surface neighbours with cluster interface atoms
        coordination = Coordination(inputfile, cluster_elements, support)
        c_coord = coordination.cluster_coordinating
        interface_distances = coordination.cluster_support_distances
        interface_indexes = coordination.support_coordinating

# gcn_i = dictionary with the generalised coordination number per atom
# c_surf = indexes of cluster atoms with coordination within the cluster lower than its bulk
        gcn = Generalised_coodination(inputfile, cluster_elements, support)
        gcn_i = gcn.gcn
        c_surf = gcn.cluster_surface_index

# c_surf_area = area exposed to the vacuum according to interpolated area/atom in the Library
        area = Areas(inputfile, cluster_elements, support)
        c_surf_area = area.cluster_surface_area

# e_coh = predicted cohesion energy in eV/atom of the whole cluster
# e_atom = sum of all the atomic energies in eV in the cluster
# e_adh = adhesion energy in eV as a function of the distance to the sites and the coordination within the cluster
# e_total = total energy in eV
# e_binding = binding energy in eV
# e_cluster_surface = predicted surface energy on the cluster atoms exposed to the vacuum.

        self.e_coh, e_atom = self.e_cohesion(system, c_coord, gcn_i)
        self.e_adh = self.e_adhesion(interface_distances, system, support, c_coord, interface_indexes)
        self.e_total = self.e_adh + e_slab + e_atom + len(c_coord) * self.e_coh
        self.e_binding = (self.e_total - e_slab - e_atom)/len(c_coord)
        self.e_cluster_surface = surface_energy(system, c_coord, c_surf, c_surf_area)

    def e_cohesion(self, system, c_coord, gcn_i):
#        cluster_elements = []
        e_coh = 0
        e_atom = 0
        average_coordination = 0
        for i in [n for n in c_coord if len(c_coord[n]) > 0]:
            average_distance = 0
            average_distance_vector = np.zeros(3)
            for j in c_coord[i]:
                average_distance += system.get_distance(int(i), int(j), mic=True)/len(c_coord[i])
                average_distance_vector += np.array(system.get_distance(int(i), int(j), mic=True, vector=True))
#                           element, cc, distance, distance, vector, gcn
            e_coh += ecoh_trend([system[int(i)].symbol], len(c_coord[i]),
                                  average_distance, list(average_distance_vector), gcn_i[int(i)])/len(c_coord)
            e_atom += float(isolated_atoms(system[int(i)].symbol))
            average_coordination += len(c_coord[i]) / len(c_coord)
#            if system[int(i)].symbol not in cluster_elements:
#                cluster_elements.append(system[int(i)].symbol)
#        e_coh = ecoh_trend(cluster_elements, average_coordination)

        return e_coh, e_atom

    def e_adhesion(self, interface_distances, system, support, c_coord, interface_indexes):
        interface_adh_e = []
        for i in interface_distances:
            v_x = [0, 0, 0]
            v_y = [0, 0, 0]
            for j in interface_indexes[i]:
                if interface_distances[i][0] == system.get_distance(int(i), int(j), mic=True):
                    v_x = system.get_distance(int(i), int(j), mic=True, vector=True)
                if interface_distances[i][1] == system.get_distance(int(i), int(j), mic=True):
                    v_y = system.get_distance(int(i), int(j), mic=True, vector=True)
#                                           support, element, icc, distance_a, distance_b, vector_distance_a, vector_distance_b
            adh_e, reference_e, e_min, distances_opt = e_adh_energies(support,
                                                                      system[int(i)].symbol,
                                                                      len(c_coord[str(i)]),
                                                                      interface_distances[i][0],
                                                                      interface_distances[i][1], v_x, v_y)
            interface_adh_e.append([i, round(adh_e, 5), round(e_min, 5),
                                    round(interface_distances[i][0]/distances_opt[0], 3),
                                    round(interface_distances[i][1]/distances_opt[1], 3)])
# Adhesion Energy Prediction RULES
# there is the distinction between two adsorption sites, i.e., strong and weak.
# interaction with the stronger site, i.e., sites[0], has preference over sites[1]
        interface_adh_e.sort(key=lambda x: x[1])
        interface_indexes = [interface_adh_e[i][0] for i in range(len(interface_adh_e))]
        primary_sites = [interface_adh_e[0][0]]
        secondary_sites = []
        for n in range(1, len(interface_adh_e)):
            i = interface_adh_e[n][0]
            if i not in secondary_sites or interface_adh_e[n][1] < interface_adh_e[n][2]*0.60: 				            # 0.60 << arbitrary parameter
                if interface_adh_e[n][1] < interface_adh_e[0][1]*0.70:												    # 0.70 << arbitrary parameter
                    primary_sites.append(i)
                    for j in c_coord[str(i)]:
                        if j in interface_indexes:
                            secondary_sites.append(j)
        if len(interface_adh_e) == len(primary_sites):
            for n in range(1, len(interface_adh_e)):
                i = interface_adh_e[n][0]
                if interface_adh_e[n][3]/interface_adh_e[n][4] > 1.0:
                    primary_sites.remove(i)
                    secondary_sites.append(i)
        secondary_sites = set([interface_adh_e[i][0] for i in range(len(interface_adh_e)) if
                               interface_adh_e[i][0] not in primary_sites])

# Predict Adhesion energy
        e_adh = 0
        for n in range(len(interface_adh_e)):
            if interface_adh_e[n][0] in primary_sites:
                e_adh += interface_adh_e[n][1]/len(primary_sites)
            elif interface_adh_e[n][0] in secondary_sites:
                e_adh += interface_adh_e[n][1]/len(secondary_sites)

        return e_adh

