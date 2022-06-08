'''

 File collecting properties data

'''

import numpy as np
from Library import sites, areas

"""
        Provides the distance from the cluster's centre of mass to the average top surface
"""


class Properties:
    def __init__(self, system, cluster_elements, support, interface_cluster_index, cluster_coordinating,
                 gcn_cluster_surface_index):

        interface_height, cluster_cm_surface_distance, cluster_mass_centre, support_cluster_min_distance = \
            self.cluster_surface_distance(system, cluster_elements, support, interface_cluster_index)
        mean_distance = self.mean_interatomic_distance(system, cluster_coordinating)
        cluster_interface_area, cluster_surface_area = self.areas(system, interface_cluster_index, cluster_coordinating,
                                                                  gcn_cluster_surface_index)
        clustering, sphericity, longest_c_cm_distance, shortest_c_cm_distance = self.sphericity(system,
                                        cluster_mass_centre, gcn_cluster_surface_index,
                                        float(cluster_interface_area + cluster_surface_area))

        values = [interface_height,               # "cs_dist" = Distance (in Å) between the surface and the cluster
                  cluster_cm_surface_distance,    # "cm_dist" = Distance (in Å) between the average surface hight and the cluster's centre of mass
                  mean_distance,                  # "cc_dist" = mean interatomic distance calculated over the first neighbour for each of the atoms in the cluster
                  longest_c_cm_distance,          # "L_c_cm" = The longest distance (in Å) between a surface atom and the cluster's centre of mass
                  shortest_c_cm_distance,         # "S_c_cm" = The shorters distance (in Å) between a surface atom and the cluster's centre of mass
                  clustering,                     # "clustering" = average difference betwenn the mean interatomic distance in the cluster and the distance between the cluster's expose and interface atoms and the centre of mass
                  sphericity,                     # "sphericity" = ratio between the cluster's area (expose and interface) and the one of a sphere with the average radius
                  cluster_interface_area,         # area coordinating the support according to interpolated area/atom in the Library
                  cluster_surface_area            # area exposed to the vacuum according to interpolated area/atom in the Library
                 ]
        keys = ["cc_dist", "cs_dist", "cm_dist", "L_radius", "S_radius", "clustering", "sphericity", "c_i_area",
                "c_s_area"]

        for i in support_cluster_min_distance:      # "dist_X" = Average of minimum distances (in Å) between the surface sites (X) and the interface clusters atom
            values.append(support_cluster_min_distance[i])
            keys.append("dist_" + i)

        self.properties = {}
        for i in range(len(keys)):
            self.properties[keys[i]] = values[i]


    def cluster_surface_distance(self, system, cluster_elements, support, c_interface):
        # c_interface = indexes of cluster atoms coordinating with a site at 1.5 * optimised distance
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]
        support_index = [system[i].index for i in range(len(system)) if system[i].symbol in sites(support)]
        support_zmax = max([system[i].position[2] for i in support_index])
        sites_index = [i.index for i in system if i.symbol in sites(support) and i.position[
            2] >= support_zmax - 1]  # gets the site atoms index in the support

        interface_height, cluster_cm_surface_distance, cluster_mass_centre = self.mass_centre(system, cluster_index,
                                                                                              sites_index, c_interface)
        support_cluster_min_distance = self.interface_distance(system, support, sites_index, c_interface)

        return interface_height, cluster_cm_surface_distance, cluster_mass_centre, support_cluster_min_distance

    def mass_centre(self, system, cluster_index, sites_index, c_interface):
        pos = [sum([(system[i].position[j] * system[i].mass) for i in cluster_index]) for j in range(3)]
        cluster_mass_centre = [float(pos[j] / sum([system[i].mass for i in cluster_index])) for j in range(3)]
        cluster_cm_surface_distance = float(cluster_mass_centre[2] - sum([system[i].position[2] for
                                                                          i in sites_index]) / len(sites_index))
        if len(c_interface) > 0:
            c_height_average = sum([system[i].position[2] for i in c_interface]) / len(c_interface)
        else:
            cluster_zmin = min([system[i].position[2] for i in cluster_index])
            c_inter_index = [i for i in cluster_index if system[i].position[2] <= cluster_zmin + 1]
            c_height_average = sum([system[i].position[2] for i in c_inter_index]) / len(c_inter_index)

        interface_height = float(
            c_height_average - sum([system[i].position[2] for i in sites_index]) / len(sites_index))

        return interface_height, cluster_cm_surface_distance, cluster_mass_centre

    def interface_distance(self, system, support, sites_index, c_interface):
        cs_distance = {}
        for site in sites(support):
            distances = []
            dist_array = []
            for n in c_interface:
                for j in sites_index:
                    dist_array.append(system.get_distance(n, j, mic=True, vector=False))
                distances.append(min(dist_array))
            cs_distance[site] = (float(sum(distances) / len(distances)))
        return cs_distance

    # Provides areas of the interface and cluster surface in A^2
    def areas(self, system, c_interface, c_coord, c_surf):
        # c_interface = indexes of cluster atoms coordinating with a site at 1.5 * optimised distance
        # c_coord = dictionary with the indexes of coordinating atoms within the cluster
        # c_surf = indexes of cluster atoms with coordination within the cluster lower than its bulk
        return self.interface(system, c_coord, c_interface), self.surface(system, c_coord, c_surf)

    def interface(self, system, c_coord, c_interface):
        return float(sum([areas(system[i].symbol, len(c_coord[str(i)])) for i in c_interface]))

    def surface(self, system, c_coord, c_surf):
        return float(sum([areas(system[i].symbol, len(c_coord[str(i)])) for i in c_surf]))

    # Provides mean interatomic distance between first neigbours in Angstroms
    def mean_interatomic_distance(self, system, c_coord):
        # c_coord = dictionary with the indexes of coordinating atoms within the cluster
        mean_distance = 0
        for i in c_coord:
            for j in c_coord[str(i)]:
                mean_distance += system.get_distance(int(i), int(j), mic=True, vector=False) / len(c_coord[str(i)])
        mean_distance = float(mean_distance / len(c_coord))
        return mean_distance

    # Provides a measure of how spherical the cluster is by comparing its expose area to the one of a perfect sphere
    def sphericity(self, system, c_mass_centre, c_surf, c_area):
        # c_mass_centre = the centre of mass from the cluster atoms considering their mass.
        # c_surf = indexes of cluster atoms with coordination within the cluster lower than its bulk
        # c_area = sum of expose and interface cluster areas
        distances = []
        for i in c_surf:
            dist_vector = np.subtract(np.array(system[i].position), c_mass_centre)
            distances.append(np.sqrt(dist_vector.dot(dist_vector)))
        sorted(distances)
        # level of sphericity
        sphericity = float(4 * np.pi * (sum(distances) / len(distances)) ** 2 / c_area)
        # level of clustering
        clustering = float((distances[-1] - distances[0]) / distances[-1])
        return clustering, sphericity, distances[-1], distances[0]
