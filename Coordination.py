"""
    by David Saadatmandi and Alberto Roldan

    Reads an input with geometries in XYZ format and calculates the
    the cluster coordination with the surface (cs)
    and coordination between metal atoms forming the cluster (cc).

"""


from ase import neighborlist
from Library import opt_atom_distance, sites, ecoh_bulk


class Coordination:
    def __init__(self, system, cluster_elements, support):
        cluster_size, cluster_ave_coordination, cluster_coordinating = self.cluster_coordination(system,
                                                                                                 cluster_elements)
        support_coordinating, sites_index_all, interface_cluster_index, interface_support_index, \
        site_cluster_coordination, cluster_support_distances = self.sites_coordination(system, cluster_elements,
                                                                                       support)
        # interface_cc_average "i_cc" = average interface atoms coordination within the cluster
        if len(interface_cluster_index) > 0:
            interface_cc_average = sum([len(cluster_coordinating[str(n)]) for n in interface_cluster_index]) /\
                                    len(interface_cluster_index)
        else:
            interface_cc_average = 0

        self.coordination = list([cluster_size,             # "N" = Total number of atoms forming the cluster
                             len(interface_cluster_index),  # "i_c" =  Number of cluster atoms at the interface
                             site_cluster_coordination,     # "cs_X" = Dictionary with the number of coordinating cluster atoms with the support sites (X)
                             interface_cc_average,          # "i_cc" = average interface atoms coordination within the cluster
                             cluster_ave_coordination       # "cc" = Average atomic coordination within the cluster
                             ])
        self.coordination_labels = list(["N", "i_c", [("cs_" + site) for site in sites(support)], "i_cc", "cc"])

        self.coordination_complete = list([cluster_size,                    # 0 "N" = Total number of atoms forming the cluster
                                            len(interface_cluster_index),   # 1 "i_c" =  Number of cluster atoms at the interface
                                            site_cluster_coordination,      # 2 "cs_X" = Dictionary with the number of coordinating cluster atoms with the support sites (X)
                                            interface_cc_average,           # 3 "i_cc" = average interface atoms coordination within the cluster
                                            cluster_ave_coordination,       # 4 "cc" = Average atomic coordination within the cluster
                                            cluster_coordinating,           # 5 dictionary with the indexes of coordinating atoms within the cluster
                                            support_coordinating,           # 6 dictionary with the indexes of coordinating sites with the cluster interface atoms
                                            sites_index_all,                # 7 dictionary with the indexes of surface sites per kind of site
                                            interface_cluster_index,        # 8 indexes of cluster atoms coordinating with a site at 1.5 * optimised distance
                                            interface_support_index,        # 9 indexes of the support coordinating with the cluster at 1.5 * optimised distance
                                            cluster_support_distances       # 10 dictionary of interface cluster atoms with the minimum distances to site X and Y.
                                        ])

    def cluster_coordination(self, system, cluster_elements):
        average_coordination = 0
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]
#        del cluster_atoms[[atom.index for atom in cluster_atoms if atom.symbol not in cluster_elements]]
        coordinating = {}
        if len(cluster_index) > 1:
            cutoff = neighborlist.natural_cutoffs(system, mult=1.3)
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
#        cs_distance = {}
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]
        cluster_zmin = min([i.position[2] for i in system if i.symbol in cluster_elements])
        interface_c_index = [i.index for i in system if i.symbol in cluster_elements and i.position[2] <= cluster_zmin + 1]
        cluster_support_distances = {}
        for site in sites(support):
#            cs_distance[site] = 0
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
            site_cluster_coordination[site] = int(len(coordinating))

        return coordinating, s_sites, interface_cluster_index, interface_support_index, \
               site_cluster_coordination, cluster_support_distances


"""
        Provides information related to the Generalised Coordination Number
        original reference: doi/10.1002/anie.201402958
"""
class Generalised_coodination:
    def __init__(self, system, cluster_elements, support):
        # c_coord = dictionary with the indexes of coordinating atoms within the cluster
        c_coord = Coordination(system, cluster_elements, support).coordination_complete[5]   # cluster_coordinating

        gcn_average, gcn, cluster_surface_index = self.cluster_gcn(system, cluster_elements, c_coord)

        self.generalised = list([gcn_average,                  # 0 "GCN" = average generalised coordination number among the atoms in the cluster
                                gcn,                           # 1 dictionary with the generalised coordination number for each atom in the cluster
                                cluster_surface_index          # 2 indexes of cluster atoms with coordination within the cluster lower than its bulk
                                ])
        self.gcn_labels = list(["GCN"])

    def cluster_gcn(self, system, cluster_elements, c_coord):
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]
        gcn_average = 0
        gcn = {}
        cluster_bulk_index = []
        cluster_surface_index = []
        for n in cluster_index:
            i_gcn = 0
            e_coh, bulk_coordination = ecoh_bulk([system[n].symbol])
            if len(c_coord[str(n)]) > 0:
                for j in c_coord[str(n)]:
                    i_gcn += len(c_coord[str(j)])/ecoh_bulk([system[j].symbol])[1]           # coordination of the coordinating atom of n
                gcn[n] = i_gcn
            if len(c_coord[str(n)]) < bulk_coordination and n not in cluster_surface_index:  # exposed atoms at the surface
                cluster_surface_index.append(n)
            if len(c_coord[str(n)]) == bulk_coordination and n not in cluster_bulk_index:
                cluster_bulk_index.append(n)
        if len([gcn[i] for i in gcn]) > 0:
            gcn_average = float(sum([gcn[i] for i in gcn])/len(gcn))

        return gcn_average, gcn, cluster_surface_index

