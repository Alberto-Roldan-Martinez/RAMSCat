"""
    Alberto Roldan

    Last Update: 14/03/2022
        create supports library
"""

import numpy as np
from scipy.misc import derivative

def e_support(support, size):                      # energies calculated with RPBE
    e_slab = 0
#                support, x, y, z, E(eV), path
    supports = {("MgO", 2, 2, 4, -368.4459, "/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Supports/MgO/MgO/Surface")}

    for sys in supports:
        if sys[0] == support and sys[3] == int(size[-1]):
            e_slab = sys[4]/(sys[1]*sys[2]) * (int(size[0])*int(size[1]))
    return e_slab

