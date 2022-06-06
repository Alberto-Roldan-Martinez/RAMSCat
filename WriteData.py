"""
    Versions:
        Alberto: 08/2019

    STRUCTURE:
        - check structure and energy (OUTCAR) files
        - check boundary conditions
        - Import coordination data
        - Import energy from the OUTCAR

"""

from ase.io import read


class Write_labels:
    def __init__(self, outfile, data):
        output = open(outfile, "w")
        output.write("#\n#\n#\n")
        for n, value in enumerate(data):
            output.write("#    Column {0:3d} = {1:^10s} = {2:}\n".format(n, value, self.label_definition(value[:3])))
        output.write("#\n#\n#\n#")
        self.length = []
        for i, value in enumerate(data):
            self.length.append(len(value) + 3)
            output.write(" {val:>{wid}}".format(wid=self.length[-1], val=value))
        output.write("\n")
        output.close()

    def label_definition(self, value):
        definitions = {
            "N": "Total number of atoms forming the cluster",
            "i_c": "Number of surface sites coordinating with the cluster",
            "cc": "Average atomic coordination within the cluster",
            "icc": "Average coordination of cluster atoms at the interface within the cluster only",
            "Ncs":  "Number of coordinating cluster atoms with the support sites (X)",
            "dis": "Average of minimum distances (in Å) between the surface sites and the clusters atoms",
            "cs_": "Distance (in Å) between the surface and the cluster",
            "cm_": "Distance (in Å) between the average surface height and the cluster's centre of mass",
            "cc_": "Mean distance (in Å) between each atom in the cluster and its first neighbours",
            "L_r": "The longest distance (in Å) between a surface atom and the cluster's centre of mass",
            "S_r": "The shortest distance (in Å) between a surface atom and the cluster's centre of mass",
            "GCN": "Average generalised coordination number for the atoms in the cluster excluding the coordination with the support",
            "c_i": "Cluster interface area (in Angstrom^2) -- check Library",
            "c_s": "Area exposed by the cluster excluding the interface (in Angstrom^2 per atom)",
            "sph": "Ratio between the cluster's area (expose and interface) and the one of a sphere with the average radius",
            "clu": "Average difference between the mean interatomic distance in the cluster and the distance between the cluster's expose and interface atoms and the centre of mass",
            "Esu": "Exposed surface energy (in J/m^2) of the cluster (not interface) -- check Library",
            "Eco": "Cohesion energy per cluster atom (in eV/atom) == (Ecluster -( N * Eatom))/N",
            "Ead": "Cluster adhesion energy (in eV) == Esystem - (Esurface + Ecluster)",
            "Eb": "Binding energy per cluster atom (in eV/atom) == (Esystem -(Esurface + N * Eatom))/N",
            "Eto": "Total energy of the system (in eV)"}

        return definitions[value]


def write_results(outfile, labels, data):
    characters_length = Write_labels(outfile, labels).length
    line = []
    for value in data:
        if type(value) is list:
            [line.append(i) for i in value]
        elif type(value) is dict:
            [line.append(value[i]) for i in value]
        else:
            line.append(value)

    output = open(outfile, "w+")
    output.write(" ")
    for n, dat in enumerate(line):
        if type(dat) is int:
            output.write(" {val:>{wid}d}".format(wid=characters_length[n], val=dat))
        elif type(dat) is str:
            output.write(" # {:<s}".format(dat))
        else:
            output.write(" {val:>{wid}.3f}".format(wid=characters_length[n], val=dat))

    output.write("\n")
    output.close()


def write_out(structure_file, energy, sphericity):
    system = read(structure_file)
    output = open("RAMSCat.out", "w+")
    output.write(" energy = {:> 12.6f} sphericity = {:>3.3f}\n"
                 .format(energy, sphericity))
    output.write(" POSITION\n---------------------------\n")  # to adapt the reading from GA

    xyz = system.get_positions()
    for i in range(len(xyz)):
        x, y, z = xyz[i]
        output.write(" {:> 5.8f}  {:> 5.8f}  {:> 5.8f}".format(x, y, z))
        output.write("\n")
    output.write("---------------------------\n")
    output.close()
