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
        coordination = Coordination(system, cluster_elements, support).coordination
        generalised = Generalised_coodination(system, cluster_elements, coordination["Others"][0]).generalised
        properties = Properties(system, cluster_elements, support, coordination["Others"][3], coordination["Others"][0],
                        generalised["Others"][1]).properties
        self.arg = [system, support, support_size, coordination["Others"][0], coordination["Others"][5],
                             coordination["Others"][1], generalised["Others"][0], generalised["Others"][1],
                             properties["c_s_area"]]

        Calculator.__init__(self)

    def calculate(self, atoms=None, properties=implemented_properties, system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        self.results = energies(*self.arg).results
