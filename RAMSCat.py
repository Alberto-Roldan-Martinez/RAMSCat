"""

This module defines an ASE interface to RAMSCat.

"""

import numpy as np
from ase.calculators.calculator import (Calculator, PropertyNotImplementedError)
from Energies import Energy_prediction as energies


class RAMSCat(Calculator):
    nolabel = True
    implemented_properties = ['energy', 'forces']
    all_changes = ['positions']

    def __init__(self, system, cluster_elements, support, support_size):
        self.arg = [system, cluster_elements, support, support_size]
        Calculator.__init__(self)

    def calculate(self, atoms=None, properties=implemented_properties, system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        self.results = energies(*self.arg).results

