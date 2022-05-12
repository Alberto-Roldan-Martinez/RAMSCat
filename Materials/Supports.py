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
    supports = {("MgO", 2, 2, 4, -368.4459, "/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Supports/MgO/MgO/Surface/2x2x4"),
#               ("MgO", 4, 4, 4, -1473.7908, "/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Supports/MgO/MgO/Surface/4x4x4"),
                ("MgO", 8, 8, 2, -715.3181, "/home/alberto/RESEARCH/OTHER/DATASET/RPBE/Supports/MgO/MgO/Surface/8x8x2"),
                }

    for sys in supports:
        if sys[0] == support and sys[1] == int(size[0]) and sys[2] == int(size[1]) and sys[3] == int(size[2]):
            e_slab = sys[4]
        else:
            print(" ERROR : The support is not defined in the library")
            exit(0)

    return e_slab

