'''
    by David Saadatmandi and Alberto Roldan

    Reads an input with geometries in XYZ format and calculates the
    the cluster coordination with the surface (cs)
    and coordination between metal atoms forming the cluster (cc).

'''


import numpy as np
from ase import neighborlist
from ase.io import read
from ase.build import bulk, molecule
from Library import opt_atom_distance, sites


class Coordination:
    def __init__(self, inputfile, cluster_elements, support):
        system = read(inputfile)

# cluster_size "N" = Total number of atoms forming the cluster
# cluster_ave_coordination "cc" = Average atomic coordination within the cluster
# cluster_coordinating = dictionary with the indexes of coordinating atoms within the cluster
        self.cluster_size, self.cluster_ave_coordination, self.cluster_coordinating = \
            self.cluster_coordination(system, cluster_elements)
        self.cluster_coord_labels = ["cc"]

# support_cluster_min_distance "dist_X" = Average of minimum distances (in Å) between the surface sites (X) and the clusters atom
# cluster_support_distances = dictionary of interface cluster atoms with the minimum distances to site X and Y.
# support_coordinating = dictionary with the indexes of coordinating sites with the cluster interface atoms
# sites_index_all = dictionary with the indexes of surface sites per kind of site
# interface_cluster_index = indexes of cluster atoms coordinating with a site at 1.5 * optimised distance
# interface_support_index = indexes of the support coordinating with the cluster at 1.5 * optimised distance
# site_cluster_coordination "cs_X" = dictionary with the number of coordinating cluster atoms with the support sites (X)
        self.support_cluster_min_distance, self.support_coordinating, self.sites_index_all, \
        self.interface_cluster_index, self.interface_support_index, self.site_cluster_coordination,\
        self.cluster_support_distances = self.sites_coordination(system, cluster_elements, support)

        self.interface_cluster = len(self.interface_cluster_index)
        self.site_cluster_coordination_label = [("cs_" + site) for site in sites(support)]
        self.support_cluster_min_distance_labels = [("dist_" + site) for site in sites(support)]

# interface_cc_average "i_cc" = average interface atoms coordination within the cluster
        if len(self.interface_cluster_index) > 0:
            self.interface_cc_average = sum([len(self.cluster_coordinating[str(n)]) for n in self.interface_cluster_index]) /\
                                    len(self.interface_cluster_index)
        else:
            self.interface_cc_average = 0

    def cluster_coordination(self, system, cluster_elements):
        average_coordination = 0
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]
#        del cluster_atoms[[atom.index for atom in cluster_atoms if atom.symbol not in cluster_elements]]
        coordinating = {}
        if len(cluster_index) > 1:
            cutoff = sum([self.bulk_distances(system[i].symbol)*1.3 for i in cluster_index])\
                     /len(cluster_index)
#            cutoff = neighborlist.natural_cutoffs(system, mult=1.3)
            a, b = neighborlist.neighbor_list('ij', system, cutoff)
            for i in cluster_index:
                coordinating[str(i)] = [b[n] for n in range(len(a)) if a[n] == i and b[n] in cluster_index]
            average_coordination = sum([len(coordinating[str(i)]) for i in cluster_index])/len(cluster_index)
        else:
            coordinating[str(cluster_index[0])] = []

        return int(len(cluster_index)), float(average_coordination), coordinating


    def sites_coordination(self, system, cluster_elements, support):
        interface_cluster_index = []
        interface_support_index = []
        site_cluster_coordination = {}
        s_sites = {}
        coordinating = {}
        cs_distance = {}
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]
        cluster_zmin = min([i.position[2] for i in system if i.symbol in cluster_elements])
        interface_c_index = [i.index for i in system if i.symbol in cluster_elements and i.position[2] <= cluster_zmin + 1]
        cluster_support_distances = {}
        for site in sites(support):
            cs_distance[site] = 0
            distances = []
            dist_array = []
            site_cluster_coordination[site] = 0
            support_zmax = max([i.position[2] for i in system if i.symbol == site])
            sites_index = [i.index for i in system if i.symbol == site and i.position[2] >= support_zmax - 1]  # gets the site atoms index in the support
            s_sites.update({site: sites_index})
            cluster_zmin = min([i.position[2] for i in system if i.symbol in cluster_elements])
            interface_cluster_index = [i.index for i in system if i.symbol in cluster_elements and i.position[2] <= cluster_zmin + 1]


            optimised_distance = [opt_atom_distance(support, site, i) for i in cluster_elements]
            distance_cutoff = sum(optimised_distance) / len(optimised_distance)
            cutoff = neighborlist.natural_cutoffs(system, distance_cutoff)
            a, b, d = neighborlist.neighbor_list('ijd', system, cutoff)

            for n in interface_c_index:
                coord = [b[i] for i in range(len(a)) if a[i] == n and b[i] in sites_index and d[i] <= distance_cutoff*1.5]
                if n not in coordinating:
                    coordinating[n] = coord
                else:
                    for i in coord:
                        coordinating[n].append(i)
                if len(coord) > 0:
                    if n not in interface_cluster_index:
                        interface_cluster_index.append(n)
                    for j in coord:
                        if j not in interface_support_index:
                            interface_support_index.append(j)
                    distances.append(min([d[i] for i in range(len(a)) if a[i] == n and b[i] in coord]))
                else:
                    for j in sites_index:
                        dist_array.append(system.get_distance(n, j, mic=True, vector=False))
                    distances.append(min(dist_array))
                if n not in cluster_support_distances:
                    cluster_support_distances[n] = [min(distances)]
                else:
                    cluster_support_distances[n].append(min(distances))
            if len(distances) > 1:
                cs_distance[site] = float(sum(distances)/len(distances))
            else:
                cs_distance[site] = float(min(distances))
            site_cluster_coordination[site] = int(len(coordinating))

        return cs_distance, coordinating, s_sites, interface_cluster_index, interface_support_index, \
               site_cluster_coordination, cluster_support_distances

    def bulk_distances(self, element):
        try:
            bulkdistance = sum(bulk(element).get_cell_lengths_and_angles()[
                               0:3]) / 3  # the average distance between atoms in the bulk
        except ValueError:
            molec = str(element.symbol + '2')
            bulkdistance = molecule(molec).get_distance(0, 1)  # interatomic distance in a diatomic molecule
        except ValueError:
            bulkdistance = 1  # minimum distance by default
        return bulkdistance


