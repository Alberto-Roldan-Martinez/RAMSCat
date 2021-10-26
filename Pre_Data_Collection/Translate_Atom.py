'''
	translates an atom long the Z axis

	USAGE: ~.py structure_file atom_index Bottom/Top (default: top)

'''

import os
import sys
from ase.io import read, write


displacements = [-0.2, -0.1, -0.05, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5, 2, 3]		# list of displacements along Z
structure_file = sys.argv[1] 		# file with the structure
atom_index = int(sys.argv[2])			# atom index to translate
try:
	atom_position = sys.argv[3]		# atom index to translate
	print(" The atom {} is sitting at the {} of the structure" .format(atom_index, atom_position))
except:
	atom_position = "Top"			# original position of atom to translate
	print(" The atom {} is sitting at the {} of the structure" .format(atom_index, atom_position))
	pass


def Rewrite(in_file_name, atom_index, displacement):
	f = open(in_file_name)
	structure = f.readlines()
	f.close()
	out_file_name = str(in_file_name[:-5]) + "_i" + str(atom_index) + "_d" + str(displacement) + ".vasp"
	out_file = open(out_file_name, 'w+')
	out_file.write("The atom {} has been translated by {} Angstroms along the Z axis\n" .format(atom_index, displacement))

	for i in range(1, 9):
		line = structure[i].split()
		if line[0].startswith("C") or line[0].startswith("D"):
			xyz_line = i + 1
	for line in range(1, xyz_line):
		out_file.write(structure[line])
		elements = [int(i) for i in structure[6].split()]
	for n_atoms in elements:
		xyz = []
		i = xyz_line				# lines in the POSCAR before the xyz positions
		for line in range(i, n_atoms+i):
			xyz_line = [float(n) for n in structure[line].split()[:3]]
			for n in structure[line].split()[3:]:
				xyz_line.append(n)
			xyz.append(xyz_line)
		i += n_atoms
		xyz = sorted(xyz, key=lambda x: x[2], reverse=False)
		for a in xyz:
			if len(a) == 3:
				out_file.write(" {:>15.11f} {:>15.11f} {:>15.11f}\n" .format(float(a[0]), float(a[1]),
																					  float(a[2])))
			else:
				out_file.write(" {:>15.11f} {:>15.11f} {:>15.11f}  {:s} {:s} {:s}\n" .format(float(a[0]), float(a[1]),
																					  float(a[2]), a[3], a[4], a[5]))
	out_file.close()
# ----------------------------------------------------------------------------------------------------------------------


for d in displacements:
	structure = read(structure_file)
	if atom_position.startswith("T") or atom_position.startswith("t"):
		structure.positions[atom_index] = structure.positions[atom_index]+[0, 0, d]
	else:
		structure.positions[atom_index] = structure.positions[atom_index]-[0, 0, d]

	write(structure_file+".vasp", structure, direct=False, vasp5=True, sort=True, ignore_constraints=False)
	Rewrite(structure_file+".vasp", atom_index, d)
	os.remove(structure_file+".vasp")

# for i in -0.2 -0.1 -0.05 0.05 0.1 0.2 0.3 0.4 0.5 1 1.5 2 3 ; do a=$(echo $i m$i |awk '{if ($1 < 0) print $2; else print $1}'); rm -rf $a; mkdir $a; cp INCAR KPOINTS run.sh $a; mv CONTCAR_i*_d$i\.vasp $a/POSCAR; cd $a; cp POSCAR POSCAR_0; sbatch run.sh; cd ..;  done;








