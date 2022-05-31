

from Library import sites
from Coordination import Coordination, Generalised_coodination
from Library import areas


"""
        Provides the distance from the cluster's centre of mass to the average top surface
"""
class Cluster_surface_distance:
    def __init__(self, system, cluster_elements, support):
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]
        support_index = [system[i].index for i in range(len(system)) if system[i].symbol in sites(support)]

        support_zmax = max([system[i].position[2] for i in support_index])
        sites_index = [i.index for i in system if i.symbol in sites(support) and i.position[2] >= support_zmax - 1]  # gets the site atoms index in the support
        pos = sum([(system[i].position[2] * system[i].mass) for i in cluster_index])
        self.cluster_mass_centre = float(pos/sum([system[i].mass for i in cluster_index]))

        self.cluster_cm_surface_distance = float(self.cluster_mass_centre - sum([system[i].position[2] for i in sites_index])\
                                           /len(sites_index))
        self.zlabels = str("Zdist")

        coordination = Coordination(system, cluster_elements, support)
        c_interface = coordination.interface_cluster_index
        if len(c_interface) > 0:
            c_height_average = sum([system[i].position[2] for i in c_interface])/len(c_interface)
        else:
            cluster_zmin = min([system[i].position[2] for i in cluster_index])
            c_inter_index = [i for i in cluster_index if system[i].position[2] <= cluster_zmin + 1]
            c_height_average = sum([system[i].position[2] for i in c_inter_index])/len(c_inter_index)

        self.interface_height = float(c_height_average - sum([system[i].position[2] for i in sites_index])/len(sites_index))


"""
        Provides areas of the interface and cluster surface in A^2
"""
class Areas:
    def __init__(self, system, cluster_elements, support):
# c_interface = indexes of cluster atoms coordinating with a site at 1.5 * optimised distance
# c_coord = dictionary with the indexes of coordinating atoms within the cluster
        coordination = Coordination(system, cluster_elements, support)
        c_interface = coordination.interface_cluster_index
        c_coord = coordination.cluster_coordinating
# c_surf = indexes of cluster atoms with coordination within the cluster lower than its bulk
        gcn = Generalised_coodination(system, cluster_elements, support)
        c_surf = gcn.cluster_surface_index

# cluster_interface_area = area coordinating the support according to interpolated area/atom in the Library
# cluster_surface_area = area exposed to the vacuum according to interpolated area/atom in the Library
        self.cluster_interface_area = self.interface(system, c_coord, c_interface)
        self.cluster_surface_area = self.surface(system, c_coord, c_surf)

    def interface(self, system, c_coord, c_interface):
        return float(sum([areas(system[i].symbol, len(c_coord[str(i)])) for i in c_interface]))

    def surface(self, system, c_coord, c_surf):
        return float(sum([areas(system[i].symbol, len(c_coord[str(i)])) for i in c_surf]))



