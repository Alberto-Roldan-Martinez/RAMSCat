"""

This module defines an ASE interface to RAMSCat.

"""

from ase.calculators.calculator import (Calculator, PropertyNotImplementedError)
from Energies import Energy_prediction as energies


class RAMSCat(Calculator):
    nolabel = True
    implemented_properties = ['energy', 'forces']
    all_changes = ['positions']

    def __init__(self, system, support, support_size, c_coord, interface_distances, interface_indexes,
                 gcn_i, c_surf, c_surf_area):
        self.arg = [system, support, support_size, c_coord, interface_distances, interface_indexes,
                    gcn_i, c_surf, c_surf_area]
        Calculator.__init__(self)

    def calculate(self, atoms=None, properties=implemented_properties, system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        self.results = energies(*self.arg).results
