'''
        provides information related to the Generalised Coordination Number
'''


from ase.io import read
from Coordination import Coordination
from Library import ecoh_bulk



class Generalised_coodination:
    def __init__(self, inputfile, cluster_elements, support):
        system = read(inputfile)
# c_coord = dictionary with the indexes of coordinating atoms within the cluster
        coordination = Coordination(inputfile, cluster_elements, support)
        c_coord = coordination.cluster_coordinating

# gcn_average = average generalised coordination number among the atoms in the cluster
# gcn = dictionary with the generalised coordination number for each atom in the cluster
# cluster_bulk_index = indexes of cluster atoms with bulk coordination
# cluster_surface_index = indexes of cluster atoms with coordination within the cluster lower than its bulk
        self.gcn_average, self.gcn, self.cluster_bulk_index, self.cluster_surface_index = \
            self.cluster_gcn(system, cluster_elements, c_coord)

    def cluster_gcn(self, system, cluster_elements, c_coord):
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]
        gcn_average = 0
        gcn = {}
        cluster_bulk_index = []
        cluster_surface_index = []
        for n in cluster_index:
            i_gcn = 0
            e_coh, bulk_coordination = ecoh_bulk(system[n].symbol)
            if len(c_coord[str(n)]) > 0:
                for j in c_coord[str(n)]:
                    i_gcn += len(c_coord[str(j)])            # coordination of the coordinating atom of n
                gcn[n] = i_gcn/bulk_coordination
            if len(c_coord[str(n)]) < bulk_coordination and n not in cluster_surface_index:  # exposed atoms at the surface
                cluster_surface_index.append(n)
            if len(c_coord[str(n)]) == bulk_coordination and n not in cluster_bulk_index:
                cluster_bulk_index.append(n)

        if len([gcn[i] for i in gcn]) > 0:
            gcn_average = float(sum([gcn[i] for i in gcn])/len(gcn))

        return gcn_average, gcn, cluster_bulk_index, cluster_surface_index