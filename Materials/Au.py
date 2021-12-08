"""
    Alberto Roldan

    Last Update: 30/11/2021
        create Au library
"""

import numpy as np


def e_isolated_atom():                      # energies of isolated atoms in vaccuo (RPBE)
    return -0.1921


def e_coh_bulk():                            # cohesion energies at bulk coordination (RPBE)
    return [-3.5650, 12]


# TO BE UPDATED!!!!!!!!!!!!
#def e_coh_trend(cc):                         # cohesion energies trend parameter (a in logarithmic equations)
#    a = 14.25097
#    b_ecoh, b_coord = e_coh_bulk()
#    return (np.log(a)/np.log(a/(a + b_coord)) - (1/np.log(a/(a + b_coord))*np.log(a+cc))) * b_ecoh
def e_coh_trend(cc, gcn):
    parameters = []
    popts = {# coordination,
        (2, ),
        (3, ),
        (4, ),
        (5, ),
        (6, ),
        (7, ),
        (8, ),
        (9, ),
        (10, )
            }
    for i, sys in enumerate(popts):
        if sys[0] == cc:
            parameters = sys[2:]
    return parameters





def surf_areas(coordination):                # Atomic areas trend using a ??? LORENTZIAN: a + b / (4 * ( x - c )**2 - d**2)
    a, b, c, d, r2 = [17.4706, 91.61526, 1.94987, 16.54962, 0.94]   # a,b,c,d and R^2 of interpolation
    return a - (b + int(coordination)**c)/d  #a + b / (4 * (int(coordination) - c)**2 - d**2)


def surf_energies(coordination):            # Atomic Surface energies interpolated using a straight line a*x + b
    a, b, r2 = [-0.00637, 0.10435, 0.71]
    return a * int(coordination) + b


def opt_support_atom_distance(support, site):          # Defines the optimised ONE atom-support distance in Ansgtroms
    optimum_distance = 0
    opt_distances = {
                    ('MgO', 'O',  2.3179),
                    ('MgO', 'Mg', 2.6857),
                    }
    for i, sys in enumerate(opt_distances):
        if sys[0] == support:
            if sys[1] == site:
                optimum_distance = sys[2]
    return optimum_distance


def e_adh_energies(support, icc):        # popt and reference_e using Trend_AdhEnergy_Sites
    parameters = []
    popts = {# support, metal, n-1	a1, a2, a_d_eq, a_r_eq, b1, b2, b_d_eq, b_r_eq, e_reference, e_min
        ('MgO',  0, 1.24123, 1.64706, 0.59456, 2.97495, 1.67188, 1.91902,  8.47915, 1.82317, -0.1787, -1.27773),
        ('MgO',  1, 1.92146, 1.72757, 1.21384, 2.18831, 2.08383, 2.17571,  0.84201, 2.25829, -0.6396, -1.35581),
        ('MgO',  2, 1.18296, 1.57988, 0.94593, 2.88265, 1.42953, 1.96775, 13.43070, 1.60617, -0.3442, -1.65847),
        ('MgO',  3, 1.23430, 1.69203, 1.04552, 2.86776, 1.50758, 2.12566, 19.78013, 1.54962, -0.5365, -1.86308),
        ('MgO',  4, 1.32852, 1.84366, 1.05141, 2.83982, 1.41161, 2.10295, 27.30569, 1.43322, -0.5228, -1.96177),
        ('MgO',  5, 1.36965, 1.80618, 1.18215, 2.65505, 1.70042, 2.12625, 12.25984, 1.63363,  0.2132, -1.91265),
        ('MgO',  6, 1.27389, 1.76288, 0.87434, 2.84338, 1.52271, 2.09858, 17.30035, 1.53972, -0.8052, -1.71807),
        ('MgO',  7, 1.25468, 1.77093, 0.86635, 2.89051, 1.45040, 2.10521, 24.11202, 1.44576, -0.8657, -1.90081),
        ('MgO',  8, 1.47385, 1.89203, 0.88521, 2.61610, 2.08310, 2.15869,  7.15395, 1.77023, -1.1498, -1.37068),
        ('MgO',  9, 1.36951, 1.78969, 1.37531, 2.56345, 2.69254, 1.59514,  5.53220, 1.83200, -1.0898, -2.83989),
        ('MgO', 10, 1.27326, 1.80154, 1.17582, 2.85178, 1.78804, 2.07628, 23.30049, 1.51652, -1.2182, -2.69864)
            }
    for i, sys in enumerate(popts):
        if sys[0] == support:
            if sys[2] == icc:
                parameters = sys[2:]
    return parameters


