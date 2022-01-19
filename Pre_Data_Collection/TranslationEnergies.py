'''
	Finds the cohesion energy of metallic systems

	USAGE: ~.py path_for_the_DFT_output_files_contains_the_energy_and_geometry index_atom_of_interest*
	 	*to get the GCN and average distance from neighbours

'''

import os, sys
from ase.io import read
from ase import neighborlist
from Library import ecoh_bulk


inputfiles = ["OUTCAR", "CONTCAR"]
path0 = os.getcwd()
try:
	source_path = "./" + sys.argv[1]				# directory containing the source information
	os.chdir(source_path)
	path = os.getcwd()
	path_name = path.split("/")[-4]+"/"+path.split("/")[-3]+"/"+path.split("/")[-2]+"/"+path.split("/")[-1]
except:
	print(" No pathDIR provided")
	exit()
try:
	i_atom = int(sys.argv[2])					# atom of interest
except:
	print(" No atom of interest provided")
	exit()


# Reading the system for its ENERGY
atoms = read(inputfiles[0], index=-1)
atoms_index = [atoms[i].index for i in range(len(atoms))]
e_atoms = atoms.get_total_energy()                      # The total energy
elements = atoms.get_chemical_symbols()                 # The elements forming it

# Reading the system for its average COORDINATION
atoms = read(inputfiles[1])
coordinating = {}
if len(atoms) > 1:
	cutoff = neighborlist.natural_cutoffs(atoms, mult=1.25)
	a, b, d = neighborlist.neighbor_list('ijd', atoms, cutoff)
	for i in atoms_index:
		coordinating[str(i)] = [b[n] for n in range(len(a)) if a[n] == i]
else:
	coordinating[str(atoms_index[0])] = 0

# Generalised coordination of the atom of interest
i_gcn = 0
if len(coordinating[str(i_atom)]) > 0:
	for j in coordinating[str(i_atom)]:
		i_gcn += len(coordinating[str(j)])/ecoh_bulk([atoms[j].symbol])[1]            # coordination of the coordinating atom to the one of interest

# the average distance between the atom of interest and first neighbours
if len(coordinating[str(i_atom)]) > 0:
	distance = sum([d[n] for n in range(len(a)) if a[n] == i_atom])/len(coordinating[str(i_atom)])
	print(i_atom, coordinating[str(i_atom)], len(coordinating[str(i_atom)]))
else:
	distance = 0

# Shortest distance to the closest atoms in the system
#distances = []
#for atom in atoms:
#	distances.append(atoms.get_distance(i_atom, atom.index, mic=True, vector=False))
#distances = sorted(distances)
#i_distance = distances[1]

# XYZ position of the atom of interest
x, y, z = atoms[i_atom].position

ifile = open("Trend_Translation.dat", 'w+')
ifile.write("# atom  coord\tgcn\tdistance(A)\tX\t Y\t  Z\tE(eV)\t\tElements\tPath\n")
ifile.write("{:5d}\t{:>3d}\t{:>3.4f}\t{:>3.4f}\t\t{:>5.4f} {:>5.4f} {:>5.4f}\t{:>3.5f}"
			.format(i_atom, len(coordinating[str(i_atom)]), i_gcn, distance, x, y, z, e_atoms))


#for i in range((len(atoms_index))):
#	ifile.write(" {:>3d}" .format(len(coordinating[str(i)])))
ifile.write("\t# {}\t{}\n" .format(set(elements), path_name))
ifile.close()

os.chdir(path0)
