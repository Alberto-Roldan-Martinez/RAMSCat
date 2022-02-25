"""
    Alberto Roldan

    Last Update: 30/11/2021
        create Au library
"""

import numpy as np
from scipy.misc import derivative

def e_isolated_atom():                      # energies of isolated atoms in vaccuo (RPBE)
    return -0.1921


def e_coh_bulk():                            # cohesion energies at bulk coordination (RPBE)
    return [-3.5650, 12]


def e_coh_trend(cc, distance, vector_distance, gcn):
    def morse_3D(x, a1, a2, a3, a_d_eq, a_r_eq, b1, b2, b3, b_d_eq, b_r_eq):
        return a_d_eq * (np.exp(-2 * a1 * (x[0] - a_r_eq)) - 2 * np.exp(-(a2 * (x[0] - (a_r_eq + a3 / x[1]))))) + \
               b_d_eq * (np.exp(-2 * b1 * (x[0] - b_r_eq)) - 2 * np.exp(-(b2 * (x[0] - (b_r_eq + b3 / x[1])))))  # MORSE potentia
    popts = {# coordination
        (0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.),
        (1, 8.27328, 30.81257, 0.00069, 1.80138, 2.01144, 1.57267, 1.33302, 0.00328, 2.02546, 2.48098),
        (2, ),
        (3,  0.30920,  0.74111, -0.38310, 2.55694, 3.32052,  2.48875, 29.43531, -2.75769, 2.51986, 2.43047),
        (4,  2.71529,  3.88079,  0.05480, 2.43175, 2.77891,  1.14273,  1.19657, -0.39842, 2.43175, 3.22909),
        (5,  1.87023, 20.38672, -1.95517, 2.51803, 2.25568,  1.87025,  0.80393, -1.09936, 3.87401, 2.29129),
        (6, ),
        (7,  1.45992,  1.64636, -0.86425, 3.78278, 3.25924, 40.86623,  3.61860,  0.20418, 3.27852, 2.72650),
        (8, ),
        (9, 38.20345,  7.92216, -2.95950, 8.42474, 2.82729,  1.80184,  0.92470, -2.47268, 5.03576, 2.70313),
        (10, 0.84165, 11.26981, -2.98870, 8.79836, 2.92734, 16.13481,  1.11706, -1.49301, 8.79836, 2.82895),
        (11, ),
        (12, )}
    for i, sys in enumerate(popts):
        if sys[0] == cc:
            arg = (gcn, ) + tuple(sys[1:])
    coh_e = morse_3D(distance, *arg)
    coh_f = vector_distance * derivative(morse_3D, distance, args=arg, dx=1e-9)

    return coh_e, coh_f
'''
- include forces with d(Morse) -->    scipy.misc.derivative(func, x0) x0=[cc, gcn]??
    2 variables: scipy.misc.derivative for multiple argument function -->stackoverflow
        forces = [ atom1, atom2, ...]
        forces of constrained atoms = 0
- implement geometry optimisation through ASE
        combine forces from cohesion and adhesion
- amend Predicting to recheck coordination every X optimisation steps if convergence hasn't been reached. x=5?? 
'''


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
    def partial_derivative(func, var=0, point=[]):
        args = point[:]
        def wraps(x):
            args[var] = x
            return func(*args)
        return derivative(wraps, point[var], dx=1e-9)

    popts = {# support,gg coordination	a1, a2, a_d_eq, a_r_eq, b1, b2, b_d_eq, b_r_eq, e_reference, e_min
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
                arguments = sys[2:-1]
                e_min = sys[-1]
    e_reference = arguments[-1]
    width = [arguments[1], arguments[5]]
    point = [distance_a, distance_b] + arguments
    e_reference = arguments[-1]
    adh_e = morse_3D(*point)
    adh_f = vector_distance_a * partial_derivative(morse_3D, 0, point=point) +\
            vector_distance_a * partial_derivative(morse_3D, 1, point=point)

    return adh_e, adh_f, e_reference, e_min, width
'''
-   add the equation to return the energy instead of the paremeters so ALSO the forces per atom can be returned --> see E-coh
'''


