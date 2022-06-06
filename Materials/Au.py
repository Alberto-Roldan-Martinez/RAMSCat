"""
    Alberto Roldan

    Last Update: 14/03/2022
        create Au library
"""

import numpy as np
import numdifftools as nd


def e_isolated_atom():                      # energies of isolated atoms in vaccuo (RPBE)
    return -0.1921


def e_coh_bulk():                            # cohesion energies at bulk coordination (RPBE)
    return [-3.5650, 12]


def e_coh_trend(cc, distance, vector_distance, gcn):
    def generalised_morse_3D(x, y, y_max, a1, a2, a3, a4, d1, d2, r1, r2, m):  # Generalised MORSE potential: https://doi.org/10.3390/en13133323
        r_eq = r1 + r2*(y_max - y)/y_max
        w = m*(y_max - y)/y_max
        if d1 < d2:
            e_dissociation = (d1/(1 + np.exp(-a4*y + a3))) + d2 		# Sigmoidal curve <<< OK
        else:
            e_dissociation = (d1/(1 + np.exp(-a4*y + a3))) - d2		# Sigmoidal curve <<< OK
        return e_dissociation * (np.exp(-2*a1*(x - r_eq)) - 2 * np.exp(-(a2+w)*(x - r_eq)))

    popts = {(0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.),
# coordination ymax,      a1,      a2,       a3,     a4,     d1,       d2,       r1,    r2,         m
# OBTAINED with 50 points interpolated from the 2D morse equations ; the R^2 are not as good but ensures a good isosurface
        (1,  0.74997, 2.35925,  1.22004, 0.83330, 0.31205, -1.89117, -2.50214, 2.37617, -0.00001,  -0.05989),
        (2,  1.50003, 0.47231, 10.38033, 1.66670, 2.77648, -2.97953, -1.47730, 2.28122, -0.05022,  -0.20637),
        (3,  2.25000, 0.55592, 11.04648, 2.50000, 3.04630, -3.11610, -1.31409, 2.48889, -0.17633,  -8.12020),
        (4,  3.19736, 0.56182,  9.84490, 2.91670, 2.42682, -3.37348, -1.77751, 2.50383, -0.04212,   1.51983),
        (5,  3.85000, 0.63111, 11.06300, 3.15602, 1.84629, -3.20573, -2.30570, 2.57143, -0.00001,  -1.26321),
        (6,  3.67497, 0.49590,  9.26037, 4.08330, 1.85846, -3.72384, -2.11259, 2.49529, -0.04408,   2.15280),
        (7,  5.24997, 0.55942,  9.43813, 5.83330, 2.74662, -4.19281, -2.06843, 2.52905, -0.17633,  -9.60421),
        (8,  6.00003, 0.47893,  6.32275, 2.83330, 5.00000, -4.50049, -2.82652, 2.45818, -0.17633,  -4.67498),
        (9,  7.45986, 0.47279,  6.84587, 7.49999, 3.03506, -4.65934, -2.62683, 2.52205, -0.26117,  -6.11175),
        (10, 7.94997, 0.68273, 12.05387, 4.33330, 5.00000, -4.02759, -2.31701, 2.68622,  0.02000, -12.12850),
        (11, 9.15003, 1.78498, 25.98537, 4.25000, 5.00000, -3.99958, -3.10658, 2.80806, -0.13959,  -6.23279),
        (12, 10.8000, 1.41918, 21.83402, 2.41176, 0.35369, -5.10972, -3.70384, 2.78977, -0.07272,   6.71450)}

    for sys in popts:
        if sys[0] == cc:
            arg = sys[1:]
    coh_e = generalised_morse_3D(distance, gcn, *arg)
    coh_f = vector_distance * nd.Gradient(generalised_morse_3D)(distance, gcn, *arg)

    return coh_e, coh_f


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


def e_adh_energies(support, icc, distance_a, distance_b, vector_distance_a, vector_distance_b):
    def morse_3D(x, y, a1, a2, a_d_eq, a_r_eq, b1, b2, b_d_eq, b_r_eq, e_reference):
        return (a_d_eq * (np.exp(-2*a1*(x - a_r_eq)) - 2 * np.exp(-a2*(x - a_r_eq*np.sin(y/x)))) +
                b_d_eq * (np.exp(-2*b1*(y - b_r_eq)) - 2 * np.exp(-b2*(y - b_r_eq*np.sin(y/x))))) + e_reference
    popts = {# support, coordination	a1, a2, a_d_eq, a_r_eq, b1, b2, b_d_eq, b_r_eq, e_reference, e_min
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
    for sys in popts:
        if str(sys[0]) == str(support):
            if int(sys[1]) == int(icc):
                arguments = [i for i in sys[2:-1]]
                e_min = sys[-1]
    width = [arguments[1], arguments[5]]
    point = [distance_a, distance_b] + arguments
    e_reference = arguments[-1]
    adh_e = morse_3D(*point)
    adh_f = np.add(vector_distance_a, vector_distance_b) * nd.Gradient(morse_3D)(distance_a, distance_b, *arguments)

    return adh_e, adh_f, e_reference, e_min, width

