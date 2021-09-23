'''
        provides areas of the interface and cluster surface in A^2
'''


from ase.io import read
from Coordination import Coordination
from GCN import Generalised_coodination
from Library import areas


class Areas:
    def __init__(self, inputfile, cluster_elements, support):
        system = read(inputfile)
# c_interface = indexes of cluster atoms coordinating with a site at 1.5 * optimised distance
# c_coord = dictionary with the indexes of coordinating atoms within the cluster
        coordination = Coordination(inputfile, cluster_elements, support)
        c_interface = coordination.interface_cluster_index
        c_coord = coordination.cluster_coordinating
# c_surf = indexes of cluster atoms with coordination within the cluster lower than its bulk
        gcn = Generalised_coodination(inputfile, cluster_elements, support)
        c_surf = gcn.cluster_surface_index

# cluster_interface_area = area coordinating the support according to interpolated area/atom in the Library
# cluster_surface_area = area exposed to the vacuum according to interpolated area/atom in the Library
        self.cluster_interface_area = self.interface(system, c_coord, c_interface)
        self.cluster_surface_area = self.surface(system, c_coord, c_surf)

    def interface(self, system, c_coord, c_interface):
        return float(sum([areas(system[i].symbol, len(c_coord[str(i)])) for i in c_interface]))

    def surface(self, system, c_coord, c_surf):
        return float(sum([areas(system[i].symbol, len(c_coord[str(i)])) for i in c_surf]))
