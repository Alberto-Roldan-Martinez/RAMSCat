"""

This module defines an ASE interface to RAMSCat.

"""

from ase.calculators.calculator import (Calculator, PropertyNotImplementedError)
from Coordination import Coordination, Generalised_coodination
from Properties import Properties
from Energies import Energy_prediction as energies


class RAMSCat(Calculator):
    nolabel = True
    implemented_properties = ['energy', 'forces']
    all_changes = ['positions']

    def __init__(self, system, cluster_elements, support, support_size):
        self.system = system
        self.cluster_elements = cluster_elements
        self.support = support
        self.support_size = support_size

        Calculator.__init__(self)

    def calculate(self, atoms=None, properties=implemented_properties, system_changes=all_changes):
        coordination = Coordination(self.system, self.cluster_elements, self.support).coordination
        generalised = Generalised_coodination(self.system, self.cluster_elements, coordination["Others"][0]).generalised
        properties = Properties(self.system, self.cluster_elements, self.support, coordination["Others"][3], coordination["Others"][0],
                        generalised["Others"][1]).properties
        self.arg = [self.system, self.support, self.support_size, coordination["Others"][0], coordination["Others"][5],
                             coordination["Others"][1], generalised["Others"][0], generalised["Others"][1],
                             properties["c_s_area"]]

        Calculator.calculate(self, atoms, properties, system_changes)

        self.results = energies(*self.arg).results
