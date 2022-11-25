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
from Library import isolated_atoms, surf_energies, ecoh_bulk, ecoh_trend, e_adh_energies, supports


class Energies:
    def __init__(self, system, e_system, cluster_elements, cluster, support, support_size, c_coord, c_surf, c_surf_area):
        # c_coord = cluster_coordinating = dictionary with the indexes of coordinating atoms within the cluster
        # c_surf = cluster_surface_index = indexes of cluster atoms with coordination within the cluster lower than its bulk
        # c_surf_area = cluster_surface_index = area exposed to the vacuum according to interpolated area/atom in the Library

        e_system = e_system.get_total_energy()
        e_slab = supports(support, support_size)
        e_cluster = cluster.get_total_energy()
        e_atoms = sum([isolated_atoms(system[i].symbol) for i in range(len(system))
                       if system[i].symbol in cluster_elements])

        e_cluster_surface = surface_energy(system, c_coord, c_surf, c_surf_area)
        cohesion = float((e_cluster - e_atoms) / len(cluster))
        adhesion = float(e_system - (e_cluster + e_slab))
        binding = float((e_system - (e_slab + e_atoms))/len(cluster))

        values = [e_cluster_surface,    # "Esurf" = surface energy (SE) in J/m^2 of the cluster atoms expossed to vacuum according to interpolated SE and area / atom in the Library
                  cohesion,             # "Ecoh" = cohesion energy per atom in the cluster in eV/atom from the DFT data
                  adhesion,             # "Eadh" = adhesion energy in eV from the DFT calculated data
                  binding,              # "Eb" = binding energy in eV from the DFT calculated data
                  float(e_system)]      # "Etotal" = total energy in eV of the supported cluster from the DFT data
        keys = list(["Esurf", "Ecoh", "Eadh", "Eb", "Etotal"])

        self.energies = {}
        for i in range(len(keys)):
            self.energies[keys[i]] = values[i]


def surface_energy(system, c_coord, c_surf, c_surf_area):
    ev_to_joules = 1.60218E-19
    e_c_surf = sum([surf_energies(system[i].symbol, len(c_coord[str(i)])) for i in c_surf])
    return float(e_c_surf * ev_to_joules / (c_surf_area * 1e-20))                  # J/m^2


class Energy_prediction:
    def __init__(self, system, support, support_size, c_coord, cluster_support_distances, interface_indexes, gcn_i, c_surf,
                 c_surf_area):
        # c_coord = dictionary with the indexes of coordinating atoms within the cluster
        # cluster_support_distances = dictionary of all cluster atoms with the minimum distances to site X and Y. --->> Updated on 18/08/2022, previouly interface_distances and changed to all the atoms in the cluster, not only the interface ones.
        # interface_indexes = support_coordinating = dictionary with surface neighbours with cluster interface atoms
        # gcn_i = gnc = dictionary with the generalised coordination number for each atom in the cluster
        # c_surf = cluster_surface_index = indexes of cluster atoms with coordination within the cluster lower than its bulk
        # c_surf_area = cluster_surface_index = area exposed to the vacuum according to interpolated area/atom in the Library

        e_slab = supports(support, support_size)

        # f_coh = predicted force in the direction between the cluster atom in question and the averaged cluster neighbours
        # e_atom = sum of ALL the atomic energies in eV in the cluster
        # f_adh = predicted force in the direction between the cluster atom in question and the surface neighbours
        e_cluster_surface = float(surface_energy(system, c_coord, c_surf, c_surf_area))
        cohesion, f_coh, e_atom = self.e_cohesion(system, c_coord, gcn_i)
        adhesion, f_adh = self.e_adhesion(cluster_support_distances, system, support, c_coord, interface_indexes)
        e_total = float(adhesion + e_slab + e_atom + len(c_coord) * cohesion)
        binding = float((e_total - (e_slab + e_atom))/len(c_coord))

        values = [e_cluster_surface,    # "Esurf" = surface energy (SE) in J/m^2 of the cluster atoms expossed to vacuum according to interpolated SE and area / atom in the Library
                  cohesion,             # "Ecoh" = cohesion energy per atom in the cluster in eV/atom from the DFT data
                  adhesion,             # "Eadh" = adhesion energy in eV from the DFT calculated data
                  binding,              # "Eb" = binding energy in eV from the DFT calculated data
                  e_total]      # "Etotal" = total energy in eV of the supported cluster from the DFT data
        keys = list(["Esurf", "Ecoh", "Eadh", "Eb", "Etotal"])

        self.energies = {}
        for i in range(len(keys)):
            self.energies[keys[i]] = values[i]

        # forces = array of atomic forces in the 3 directions [[x1, y1, z1], [x2,...],...]
        forces = []
        for i in range(len(system)):
            f = np.zeros(3,)
            if str(system[i].index) in f_coh:
                # force control to smoothen large forces
                f_control = []
                for n in f_coh[str(i)]:
                    if np.abs(n) > 100:
                        f_control.append(n * 1E-2)
                    else:
                        f_control.append(n)
                f += [n for n in f_control]
#                f += [n for n in f_coh[str(i)]]

            if str(system[i].index) in f_adh:
                # force control to smoothen large forces
                f_control = []
                for n in f_adh[str(i)]:
                    if np.abs(n) > 100:
                        f_control.append(n * 1E-2)
                    else:
                        f_control.append(n)
                f += [n for n in f_control]
#                f += [n for n in f_adh[str(i)]]

            forces.append(f)

        self.results = {'energy': e_total,
                        'forces': np.array(forces),
                        'stress': np.zeros(6), 'dipole': np.zeros(3), 'charges': np.zeros(len(system)),
                        'magmom': 0.0, 'magmoms': np.zeros(len(system))}


    def e_cohesion(self, system, c_coord, gcn_i):
        e_cohesion = 0
        f_cohesion = {}
        e_atom = 0
        average_coordination = 0
        # Energy and Forces of each atom in cluster with coordination > 0 (i.e. with gcn) against each individual atom in the cluster
        indexes = []
        coh_e_library = {}
        for i in [n for n in c_coord if len(c_coord[n]) > 0]:       # atoms without neighbours don't have Coh
            coh_e = []
            f_cohesion[str(i)] = np.zeros(3)
            for j in [n for n in c_coord if n != i]:
                distance = float(system.get_distance(int(i), int(j), mic=True))
                distance_vector = np.array((system.get_distance(int(i), int(j), mic=True, vector=True)))
                #      element, cc, distance, distance, vector, gcn
                e, f = ecoh_trend([system[int(i)].symbol], len(c_coord[i]), distance, list(distance_vector),
                                  gcn_i[int(i)])
                coh_e.append([i, j, distance, e])
                f_cohesion[str(i)] += f/(len(c_coord) - 1)
            e_atom += float(isolated_atoms(system[int(i)].symbol))
            average_coordination += len(c_coord[i]) / len(c_coord[str(i)])

            # The energy is defined by the closest to atom i
            coh_e.sort(key=lambda x: x[3])
            indexes.append(str(coh_e[0][0]))
            indexes.append(str(coh_e[0][1]))
            coh_e_library[str(i)] = coh_e[0][3]
        for i in coh_e_library:
            e_cohesion += coh_e_library[i]/indexes.count(i)             # double counting i --> j and j --> i

        return float(e_cohesion/len(c_coord)), f_cohesion, float(e_atom)


    # Adhesion energy measured from ALL cluster atoms.
    # The Adhesion, should be measured for all the atoms in the cluster independently if they are in
    # contact with the surface or not.
    def e_adhesion(self, cluster_support_distances, system, support, c_coord, interface_indexes):
        interface_adh_e = []
        adh_f = {}
        support_height = []
        cluster_interface_height = []
        for i in interface_indexes:
            support_height.append(sum([system[j].position[2] for j in interface_indexes[i]])/len(interface_indexes[i]))
            cluster_interface_height.append(system[i].position[2])
        support_height_average = round(sum(support_height) / len(interface_indexes), 5)
        cluster_interface_height_ave = round(sum(cluster_interface_height)/len(interface_indexes) -
                                                 support_height_average, 5)
        for i in cluster_support_distances:             # ALL cluster atoms
            site_symbol_a, site_index_a, distance_a = cluster_support_distances[i][0]
            vector_distance_a = system.get_distance(int(i), int(site_index_a), mic=True, vector=True)
            site_symbol_b, site_index_b, distance_b = cluster_support_distances[i][1]
            vector_distance_b = system.get_distance(int(i), int(site_index_b), mic=True, vector=True)
            #                               support, element, icc, distance_a, distance_b, vector_distance_a, vector_distance_b
            try:
                e_adh, f_adh, reference_e, e_min, distances_opt = e_adh_energies(support,
                                                                             system[int(i)].symbol,
                                                                             len(c_coord[str(i)]),
                                                                             distance_a,
                                                                             distance_b,
                                                                             vector_distance_a,
                                                                             vector_distance_b)
            except:
                e_adh = 0
                e_min = 0
                continue

            interface_adh_e.append([i, round(e_adh, 5), round(e_min, 5),
                                    round(system[i].position[2] - support_height_average, 5)])
#                                    round(distance_a, 3),
#                                    round(distance_a/distances_opt[0], 3),
#                                    round(distance_b, 3),
#                                    round(distance_b/distances_opt[1], 3)])
            adh_f[str(i)] = f_adh

        interface_adh_e.sort(key=lambda x: x[3])
#        primary = []
        secondary = []
        primary_energy = 0
        secondary_energy = 0
        for index in interface_adh_e:
            if index[0] not in secondary and index[3] <= cluster_interface_height_ave and index[0] in interface_indexes:
#                primary.append(index[0])
                primary_energy += float(index[1])
                for i in c_coord[str(index[0])]:
                    if i in interface_indexes:
                        secondary.append(i)
                        i_neighbours = [k for k in c_coord[str(index[0])] if k in interface_indexes]
                        secondary_energy += sum([j[1] for j in interface_adh_e if j[0] == i])/len(i_neighbours)
            elif index[0] not in secondary and index[0] in interface_indexes:
#                primary.append(index[0])
                primary_energy += float(index[1]/(1 + len([i for i in c_coord[str(index[0])] if i in interface_indexes])))

        adh_e = primary_energy # - secondary_energy

        return float(adh_e), adh_f



