"""
    Alberto Roldan

    STRUCTURE:
        - Library of interatomic distances (metals)
            and bulk lattices as distance between sites (oxides)

"""

import numpy as np
import importlib




def isolated_atoms(element):                     # energies of isolated atoms in vaccuo (RPBE)
    return importlib.import_module("Materials."+element).e_isolated_atom()
'''  
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
'''


def ecoh_bulk(element):                         # cohesion energies at bulk coordination (RPBE)
    if type(element) is list and len(element) == 1:
        return importlib.import_module("Materials."+str(element[0])).e_coh_bulk()
    else:
        print(" --- Bulk alloys not implemented in the Library yet ---")
        exit()
'''  
   if type(element) is list and len(element) == 1:
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
        b_ecoh = ecoh[str(element[0])]
    else:
        print(" --- Bulk alloys not implemented in the Library yet ---")
        exit()

    return b_ecoh
'''
def ecoh_trend(element, cc, distance, vector_distance, gcn):                         # cohesion energies trend parameter (a in logarithmic equations)
    if type(element) is list and len(element) == 1:
        coh_e, coh_f = importlib.import_module("Materials."+str(element[0])).e_coh_trend(cc, distance,
                                                                                         vector_distance, gcn)
        return coh_e
    else:
        print(" --- Bulk alloys not implemented in the Library yet ---")
        exit()


def areas(element, coordination):                # Atomic areas trend using a ??? LORENTZIAN: a + b / (4 * ( x - c )**2 - d**2)
     return importlib.import_module("Materials."+element).surf_areas(coordination)
'''  
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

'''
def surf_energies(element, coordination):        # Atomic Surface energies interpolated using a straigh line a*x + b
     return importlib.import_module("Materials."+element).surf_energies(coordination)
'''  
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
'''


def sites(support):         # defines the common adsorption sites on a support sorted by adh, i.e., stronger to weaker.
    sites = {
        "MgO": ["O", "Mg"]
    }
    return sites[str(support)]


def opt_atom_distance(support, site, element):          # Defines the optimised ONE atom-support distance in Angstroms
    return importlib.import_module("Materials."+element).opt_support_atom_distance(support, site)
'''
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

    return optimum_distance
'''


def e_adh_energies(support, element, icc, x, y, v_x, v_y):
    adh_e, adh_f, e_reference, e_min, widths = importlib.import_module("Materials."+element).e_adh_energies(support,
                                                                                                            icc, x, y,
                                                                                                            v_x, v_y)
    return adh_e, e_reference, e_min, widths

def supports(support, size):
#   assuming that the energy of a surface is scalable with its size
    return importlib.import_module("Materials.Supports").e_support(support, size)
