'''
	Finds the cohesion energy of metallic systems

	USAGE: ~.py path_for_the_DFT_output_files_contains_the_energy_and_geometry
		OPTIONAL: index_atom_of_interest to get the GCN and average distance from neighbours

'''

import os, sys
from ase.io import read
from ase import neighborlist
from Library import isolated_atoms, ecoh_bulk


path0 = os.getcwd()
try:
	source_path = "./" + sys.argv[1]				# directory containing the source information
	os.chdir(source_path)
	try:
		i_atom = int(sys.argv[2])					# atom of interest
	except:
		i_atom = ""
		pass
except:
	pass

inputfiles = ["OUTCAR", "CONTCAR"]
path = os.getcwd()
path_name = path.split("/")[-4]+"/"+path.split("/")[-3]+"/"+path.split("/")[-2]+"/"+path.split("/")[-1]

# Reading the system for its ENERGY
atoms = read(inputfiles[0], index=-1)
atoms_index = [atoms[i].index for i in range(len(atoms))]
e_atoms = atoms.get_total_energy()                      # The total energy
elements = atoms.get_chemical_symbols()                 # The elements forming it
cohesion_e = (e_atoms - sum([isolated_atoms(i) for i in elements])) / len(elements)
cohesion_e_bulk = sum([ecoh_bulk([i])[0] for i in elements]) / len(elements)

# Reading the system for its average COORDINATION
atoms = read(inputfiles[1])
coordinating = {}
if len(atoms_index) > 1:
	cutoff = neighborlist.natural_cutoffs(atoms, mult=1.3)
	a, b, d = neighborlist.neighbor_list('ijd', atoms, cutoff)
	for i in atoms_index:
		coordinating[str(i)] = [b[n] for n in range(len(a)) if a[n] == i]
else:
	coordinating[str(atoms_index[0])] = 0
average_coordination = sum([len(coordinating[str(i)]) for i in atoms_index])/len(atoms_index)
coord_bulk = sum([ecoh_bulk([i])[1] for i in elements]) / len(elements)

# Generalised coordination of the atom of interest
if type(i_atom) is int:
	#and len(coordinating[str(i_atom)]) > 0:
	i_gcn = 0
	for j in coordinating[str(i_atom)]:
		i_gcn += len(coordinating[str(j)])            # coordination of the coordinating atom to the one of interest
	gcn = i_gcn/ecoh_bulk([atoms[i_atom].symbol])[1]
# Shortest distance to the closest atoms in the system
	distances = []
	for atom in atoms:
		distances.append(atoms.get_distance(i_atom, atom.index, mic=True, vector=False))
	distances = sorted(distances)
	distance = distances[1]
else:
	gcn = 0
	distance = 0

ifile = open("Trend_CohEnergy.dat", 'w+')
ifile.write("# ave_coord\tgcn\tdistance(A)\tE_Coh (eV.atom\N{SUPERSCRIPT minus}\N{SUPERSCRIPT ONE})\tE_Coh^Bulk\t\tElements\tPath\n")
ifile.write("{:>3.4f}\t\t{:>3.4f}\t{:>3.4f}\t\t{:>5.4f}\t\t\t{:>3.2f} {:>5.4f}\t" .format(average_coordination, gcn,
																						distance, cohesion_e,
																						coord_bulk, cohesion_e_bulk))
#for i in range((len(atoms_index))):
#	ifile.write(" {:>3d}" .format(len(coordinating[str(i)])))
ifile.write("\t# {}\t{}\n" .format(set(elements), path_name))
ifile.close()

os.chdir(path0)
