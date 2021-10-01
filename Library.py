"""
    Alberto Roldan

    STRUCTURE:
        - Library of interatomic distances (metals)
            and bulk lattices as distance between sites (oxides)

"""

import numpy as np


def opt_atom_distance(support, site, element):          # Defines the optimised ONE atom-support distance in Ansgtroms

# add also the adhesion energy of a single atom so to have eadh (0)

    optdistances = {
                    ('MgO', 'O',    'Co',   ),    # 3rd row
                    ('MgO', 'Mg',   'Co',   ),
                    ('MgO', 'O',    'Ni',   ),
                    ('MgO', 'Mg',   'Ni',   ),            
                    ('MgO', 'O',    'Cu',   2.0184),
                    ('MgO', 'Mg',   'Cu',   2.7037),    # 3rd row
                    ('MgO', 'O',    'Ru',   2.0032),    # 4th row
                    ('MgO', 'Mg',   'Ru',   2.6564),
                    ('MgO', 'O',    'Rh',   2.0836),
                    ('MgO', 'Mg',   'Rh',   2.6412),
                    ('MgO', 'O',    'Pd',   2.1020),
                    ('MgO', 'Mg',   'Pd',   2.5317),
                    ('MgO', 'O',    'Ag',   2.4385),
                    ('MgO', 'Mg',   'Ag',   2.8251),    # 4th row
                    ('MgO', 'O',    'Ir',   2.0184),    # 5th row
                    ('MgO', 'Mg',   'Ir',   2.8757),
                    ('MgO', 'O',    'Pt',   1.9885),
                    ('MgO', 'Mg',   'Pt',   2.5948),
                    ('MgO', 'O',    'Au',   2.3179),
                    ('MgO', 'Mg',   'Au',   2.6857),
                    ('C', 'C',   'Au',   2.6857),    # 5th row
                    }
    for i, sys in enumerate(optdistances):
        if sys[0] == support:
            if sys[1] == site:
                if sys[2] == element:
                    optimum_distance = sys[3]
#                   e_adh_0 = sys[4]
#   return optimum_distance, e_adh_0
    return optimum_distance


def isolated_atoms(element):                     # energies of isolated atoms in vaccuo (RPBE)
    elements = {
                 'Fe': -3.4949,                 # 3rd row
                 'Co': -2.1442,
                 'Ni': -0.8036,
                 'Cu': -0.0099,                 # 3rd row
                 'Ru': -2.5726,                 # 4th row
                 'Rh': -1.5157,
                 'Pd': -1.5043,
                 'Ag': -0.1921,                 # 4th row
                 'Ir': -1.2597,                 # 5th row
                 'Pt': -1.6193,
                 'Au': -0.1921,                 # 5th row
               }

    return elements[element]


def ecoh_bulk(element):                         # cohesion energies at bulk coordination (RPBE)
    ecoh = {
            'Co': [-7.1245, 12],    # hcp       # 3rd row
            'Ni': [-5.6273, 12],
            'Cu': [-4.0698, 12],                # 3rd row
            'Ru': [-9.4469, 12],    # hcp       # 4rd row
            'Rh': [-7.5247, 12],
            'Pd': [-5.5162, 12],
            'Ag': [-3.0364, 12],                # 4rd row
            'Ir': [-9.4589, 12],                # 5rd row
            'Pt': [-6.5738, 12],
            'Au': [-3.5650, 12],                # 5th row
           }
    return ecoh[element]

def ecoh_trend(element, cc):                         # cohesion energies trend parameter (a in logarithmic equations)
    coh_parameter = {
            'Au': 10.88828,                # 5th row
           }
    a = coh_parameter[element]
    return (np.log(a)/np.log(a/(a + ecoh_bulk(element)[1])) -
            (1/np.log(a/(a + ecoh_bulk(element)[1])))*np.log(a+cc)) * ecoh_bulk(element)[0]


def areas(element, coordination):                # Atomic areas previously calculated from surfaces as a function of the atom's coordination [0-->12]
    area = {                                    #   interpolated using a  LORENTZIAN function(x,a,b,c,d):   a + b / (4 * ( x - c )**2 - d**2)
#           element || areas vs coordination [0-->12] in Angstroms || a,b,c,d and R^2 of interpolation
# 3rd row           a        b       c       d       R^2
            'Co': [-7.1245, 12], # hcp
            'Ni': [-5.6273, 12],
            'Cu': [-4.0698, 12],
# 4rd row           a        b       c       d       R^2
            'Ru': [-9.4469, 12], # hcp
            'Rh': [-7.5247, 12],
            'Pd': [-5.5162, 12],
            'Ag': [-3.0364, 12],
# 5rd row           a        b       c       d       R^2
            'Ir': [-9.4589, 12],                
            'Pt': [-6.5738, 12],
            'Au': [17.4706, 91.61526, 1.94987, 16.54962, 0.94],
           }
    a, b, c, d, r2 = area[element]
    return a - (b + int(coordination)**c)/d  #a + b / (4 * (int(coordination) - c)**2 - d**2)


def surf_energies(element, coordination):        # Atomic Surface energies previously calculated from surfaces as a function of the atom's coordination [0-->12]
    surf_energy = {                                 #   interpolated using a straigh line a*x + b
#           element || energies vs coordination [0-->12] in Angstroms || a,b and R^2 of interpolation
# 3rd row           a        b       R^2
            'Co': [-7.1245, 12], # hcp
            'Ni': [-5.6273, 12],
            'Cu': [-4.0698, 12],
# 4rd row           a        b        R^2
            'Ru': [-9.4469, 12], # hcp
            'Rh': [-7.5247, 12],
            'Pd': [-5.5162, 12],
            'Ag': [-3.0364, 12],
# 5rd row           a      b      R^2
            'Ir': [-9.4589, 12],                
            'Pt': [-6.5738, 12],
            'Au': [-0.00637, 0.10435, 0.71],
            }
    a, b, r2 = surf_energy[element]
    return a * int(coordination) + b


def sites(support):         # defines the common adsorption sites on a support sorted by adh, i.e., stronger to weaker.
    sites = {
        "MgO": ["O", "Mg"]
    }
    return sites[str(support)]


def morse_3D_energies(support, element, icc, x, y):
    popts = {# support, metal, n-1											# popt and reference_e using Trend_AdhEnergy_Sites
#       add the icc = 0 sum of O-sites and Mg-site
        ('MgO', 'Au',  1, 1.92146, 1.72757, 1.21384, 2.18831, 2.08383, 2.17571,  0.84201, 2.25829, -0.6396, -1.35581),		# 2 atoms
        ('MgO', 'Au',  2, 1.18296, 1.57988, 0.94593, 2.88265, 1.42953, 1.96775, 13.43070, 1.60617, -0.3442, -1.65847),		# 3 atoms
        ('MgO', 'Au',  3, 1.23430, 1.69203, 1.04552, 2.86776, 1.50758, 2.12566, 19.78013, 1.54962, -0.5365, -1.86308),		# 4 atoms
        ('MgO', 'Au',  4, 1.32852, 1.84366, 1.05141, 2.83982, 1.41161, 2.10295, 27.30569, 1.43322, -0.5228, -1.96177),		# 5 atoms
        ('MgO', 'Au',  5, 1.36965, 1.80618, 1.18215, 2.65505, 1.70042, 2.12625, 12.25984, 1.63363,  0.2132, -1.91265),		# 6 atoms
        ('MgO', 'Au',  6, 1.27389, 1.76288, 0.87434, 2.84338, 1.52271, 2.09858, 17.30035, 1.53972, -0.8052, -1.71807),		# 7 atoms
        ('MgO', 'Au',  7, 1.25468, 1.77093, 0.86635, 2.89051, 1.45040, 2.10521, 24.11202, 1.44576, -0.8657, -1.90081),		# 8 atoms
        ('MgO', 'Au',  8, 1.47385, 1.89203, 0.88521, 2.61610, 2.08310, 2.15869,  7.15395, 1.77023, -1.1498, -1.37068),		# 9 atoms
        ('MgO', 'Au',  9, 1.36951, 1.78969, 1.37531, 2.56345, 2.69254, 1.59514,  5.53220, 1.83200, -1.0898, -2.83989),		# 10 atoms
        ('MgO', 'Au', 10, 1.27326, 1.80154, 1.17582, 2.85178, 1.78804, 2.07628, 23.30049, 1.51652, -1.2182, -2.69864)		# 11 atom
            }
    for i, sys in enumerate(popts):
        if sys[0] == support:
            if sys[1] == element:
                if sys[2] == icc:
                    support, element, icc, a1, a2, a_d_eq, a_r_eq, b1, b2, b_d_eq, b_r_eq, e_reference, e_min = sys
                    adh_e = (a_d_eq * (np.exp(-2*a1*(x - a_r_eq)) - 2 * np.exp(-a2*(x - a_r_eq*np.sin(y/x)))) +
                               b_d_eq * (np.exp(-2*b1*(y - b_r_eq)) - 2 * np.exp(-b2*(y - b_r_eq*np.sin(y/x))))) +\
                               e_reference		# MORSE potential

    return adh_e, e_reference, e_min, [a2, b2]




