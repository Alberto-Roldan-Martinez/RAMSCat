'''
	Finds the Adhesion energy of metallic systems on supports

	USAGE: ~.py path_for_the_DFT_output_files
		FILE_1: contains energy
		FILE_2: contains geometry

'''

import os
from ase.io import read
from ase import neighborlist
from Library import isolated_atoms, ecoh_bulk


inputfiles = ["OUTCAR", "CONTCAR"]
path = os.getcwd()
path_name = path.split("/")[-4]+"/"+path.split("/")[-3]+"/"+path.split("/")[-2]+"/"+path.split("/")[-1]

# Reading the system for its ENERGY
atoms = read(inputfiles[0], index=-1)
atoms_index = [atoms[i].index for i in range(len(atoms))]
e_atoms = atoms.get_total_energy()                      # The total energy
elements = atoms.get_chemical_symbols()                 # The elements forming it
cohesion_e = (e_atoms - sum([isolated_atoms(i) for i in elements])) / len(elements)
cohesion_e_bulk = sum([ecoh_bulk(i)[0] for i in elements]) / len(elements)

# Reading the system for its average COORDINATION
#atoms = read(inputfiles[1])
coordinating = {}
if len(atoms_index) > 1:
	cutoff = neighborlist.natural_cutoffs(atoms, mult=1.25)
	a, b = neighborlist.neighbor_list('ij', atoms, cutoff)
	for i in atoms_index:
		coordinating[str(i)] = [b[n] for n in range(len(a)) if a[n] == i]
else:
	coordinating[str(atoms_index[0])] = 0
average_coordination = sum([len(coordinating[str(i)]) for i in atoms_index])/len(atoms_index)
coord_bulk = sum([ecoh_bulk(i)[1] for i in elements]) / len(elements)

ifile = open("Trend_CohEnergy.dat", 'w+')
ifile.write("# ave_coord\tE_Coh (eV.atom\N{SUPERSCRIPT minus}\N{SUPERSCRIPT ONE})\tE_Coh^Bulk\t\tElements\tPath\n")
ifile.write("{:>3.4f}\t\t{:>5.4f}\t\t\t{:>3.2f} {:>5.4f}\t" .format(average_coordination, cohesion_e, coord_bulk, cohesion_e_bulk))
#for i in range((len(atoms_index))):
#	ifile.write(" {:>3d}" .format(len(coordinating[str(i)])))
ifile.write("\t# {}\t\t{}\n" .format(atoms.symbols, path_name))
ifile.close()
