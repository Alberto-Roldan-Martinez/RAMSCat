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
from Coordination import Coordination, Generalised_coodination
from Properties import Areas
from Library import isolated_atoms, surf_energies, ecoh_trend, e_adh_energies, supports


class Energies:
    def __init__(self, system, e_system, cluster_elements, cluster, support, support_size):
        e_system = e_system.get_total_energy()
        e_slab = supports(support, support_size)
        e_cluster = cluster.get_total_energy()
        e_atoms = sum([isolated_atoms(system[i].symbol) for i in range(len(system))
                       if system[i].symbol in cluster_elements])

# c_coord = dictionary with the indexes of coordinating atoms within the cluster
        coordination = Coordination(system, cluster_elements, support)
        c_coord = coordination.cluster_coordinating
# c_surf = indexes of cluster atoms with coordination within the cluster lower than its bulk
        gcn = Generalised_coodination(system, cluster_elements, support)
        c_surf = gcn.cluster_surface_index
# c_surf_area = area exposed to the vacuum according to interpolated area/atom in the Library
        area = Areas(system, cluster_elements, support)
        c_surf_area = area.cluster_surface_area

# cohesion = cohesion energy per atom in the cluster in eV/atom from the DFT data
# adhesion = adhesion energy in eV from the DFT calculated data
# binding = binding energy in eV from the DFT calculated data
# e_total = total energy in eV of the supported cluster from the DFT data
# e_cluster_surface = surface energy (SE) in J/m^2 of the cluster atoms expossed to vacuum according to interpolated SE and area / atom in the Library
        self.cohesion = float((e_cluster - e_atoms) / len(cluster))
        self.adhesion = float(e_system - (e_cluster + e_slab))
        self.binding = float((e_system - (e_slab + e_atoms))/len(cluster))
        self.e_total = float(e_system)
        self.e_cluster_surface = surface_energy(system, c_coord, c_surf, c_surf_area)


def surface_energy(system, c_coord, c_surf, c_surf_area):
    ev_to_joules = 1.60218E-19
    e_c_surf = sum([surf_energies(system[i].symbol, len(c_coord[str(i)])) for i in c_surf])
    return float(e_c_surf * ev_to_joules / (c_surf_area * 1e-20))                  # J/m^2


class Energy_prediction:
    def __init__(self, system, cluster_elements, support, support_size):
        e_slab = supports(support, support_size)
# c_coord = dictionary with the indexes of coordinating atoms within the cluster
# interface_distances = dictionary of interface cluster atoms with the minimum distances to site X and Y
# interface_indexes = dictionary with surface neighbours with cluster interface atoms
        coordination = Coordination(system, cluster_elements, support)
        c_coord = coordination.cluster_coordinating
        interface_distances = coordination.cluster_support_distances
        interface_indexes = coordination.support_coordinating

# gcn_i = dictionary with the generalised coordination number per atom
# c_surf = indexes of cluster atoms with coordination within the cluster lower than its bulk
        gcn = Generalised_coodination(system, cluster_elements, support)
        gcn_i = gcn.gcn
        c_surf = gcn.cluster_surface_index

# c_surf_area = area exposed to the vacuum according to interpolated area/atom in the Library
        area = Areas(system, cluster_elements, support)
        c_surf_area = area.cluster_surface_area

# e_coh = predicted cohesion energy in eV/atom of the whole cluster
# f_coh = predicted force in the direction between the cluster atom in question and the averaged cluster neighbours
# e_atom = sum of all the atomic energies in eV in the cluster
# e_adh = adhesion energy in eV as a function of the distance to the sites and the coordination within the cluster
# f_adh = predicted force in the direction between the cluster atom in question and the surface neighbours
# e_total = total energy in eV
# e_binding = binding energy in eV
# e_cluster_surface = predicted surface energy on the cluster atoms exposed to the vacuum.

        self.cohesion, f_coh, e_atom = self.e_cohesion(system, c_coord, gcn_i)
        self.adhesion, f_adh = self.e_adhesion(interface_distances, system, support, c_coord, interface_indexes)
        self.e_total = float(self.adhesion + e_slab + e_atom + len(c_coord) * self.cohesion)
        self.binding = float((self.e_total - (e_slab + e_atom))/len(c_coord))
        self.e_cluster_surface = float(surface_energy(system, c_coord, c_surf, c_surf_area))

# forces = array of atomic forces in the 3 directions [[x1, y1, z1], [x2,...],...]
        forces = []
        for i in range(len(system)):
            f = np.zeros(3,)
            if str(system[i].index) in f_coh:
                f += [i for i in f_coh[str(i)]]

            if str(system[i].index) in f_adh:
                f += [i for i in f_adh[str(i)]]
            forces.append(f)

        self.results = {'energy': self.e_total,
                        'forces': np.array(forces),
                        'stress': np.zeros(6), 'dipole': np.zeros(3), 'charges': np.zeros(len(system)),
                        'magmom': 0.0, 'magmoms': np.zeros(len(system))}

    def e_cohesion(self, system, c_coord, gcn_i):
#        cluster_elements = []
        e_cohesion = 0
        f_cohesion = {}
        e_atom = 0
        average_coordination = 0
        for i in [n for n in c_coord if len(c_coord[n]) > 0]:
            average_distance = 0
            average_distance_vector = np.zeros(3)
            for j in c_coord[i]:
                average_distance += float(system.get_distance(int(i), int(j), mic=True)/len(c_coord[i]))
                average_distance_vector += np.array((system.get_distance(int(i), int(j), mic=True, vector=True))
                                                    /len(c_coord[i]))
#                           element, cc, distance, distance, vector, gcn
            e, f = ecoh_trend([system[int(i)].symbol], len(c_coord[i]),
                              average_distance, list(average_distance_vector), gcn_i[int(i)])
            e_cohesion += e/len(c_coord)
            f_cohesion[str(i)] = f
            e_atom += float(isolated_atoms(system[int(i)].symbol))
            average_coordination += len(c_coord[i]) / len(c_coord)

        return float((e_cohesion - e_atom)/len(c_coord)), f_cohesion, float(e_atom)

    def e_adhesion(self, interface_distances, system, support, c_coord, interface_indexes):
        interface_adh_e = []
        adh_f = {}
        for i in interface_distances:
            v_x = [0, 0, 0]
            v_y = [0, 0, 0]
            for j in interface_indexes[i]:
                if interface_distances[i][0] == system.get_distance(int(i), int(j), mic=True):
                    v_x = system.get_distance(int(i), int(j), mic=True, vector=True)
                if interface_distances[i][1] == system.get_distance(int(i), int(j), mic=True):
                    v_y = system.get_distance(int(i), int(j), mic=True, vector=True)
#                                           support, element, icc, distance_a, distance_b, vector_distance_a, vector_distance_b
            e_adh, f_adh, reference_e, e_min, distances_opt = e_adh_energies(support,
                                                                      system[int(i)].symbol,
                                                                      len(c_coord[str(i)]),
                                                                      interface_distances[i][0],
                                                                      interface_distances[i][1], v_x, v_y)
            adh_f[str(i)] = f_adh
            interface_adh_e.append([i, round(e_adh, 5), round(e_min, 5),
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

# Predict Adhesion energy & forces
        adh_e = 0
        for n in range(len(interface_adh_e)):
            if interface_adh_e[n][0] in primary_sites:
                adh_e += interface_adh_e[n][1]/len(primary_sites)
            elif interface_adh_e[n][0] in secondary_sites:
                adh_e += interface_adh_e[n][1]/len(secondary_sites)

        return float(adh_e), adh_f


