'''

	Compares the columns in two input files if the values in column X and Y are < 0.25
	USAGE: ~.py a.dat b.dat

'''

import sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

icolour = ["b", "r", "k", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] # n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))] # n=7


def get_data(dataset):
	ifile = open(dataset).readlines()
	dataset = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	labels_line = len([1 for i in range(len(ifile)) if ifile[i].startswith("#") is True]) - 1
	un_labels = ifile[labels_line].split()
#     Column   0 =     N      = Total number of atoms forming the cluster
#     Column   1 =    i_c     = Number of cluster atoms coordinating the surface
#     Column   2 =    cs_O    = Number of surface sites coordinating with the cluster
#     Column   3 =   cs_Mg    = Number of surface sites coordinating with the cluster
#     Column   4 =    i_cc    = Average coordination of cluster atoms at the interface within the cluster only
#     Column   5 =     cc     = Average atomic coordination within the cluster
#     Column   6 =   dist_O   = Average of minimum distances (in Å) between the surface sites and the clusters atoms
#     Column   7 =  dist_Mg   = Average of minimum distances (in Å) between the surface sites and the clusters atoms
#     Column   8 = cs_height  = Distance (in Å) between the surface and the cluster
#     Column   9 =   Zdist    = Distance (in Å) between the average surface high and the cluster's centre of mass
#     Column  10 =    GCN     = Average generalised coordination number for the atoms in the cluster excluding the coordination with the support
#     Column  11 =  c_i_area  = Cluster interface area (in Angstrom^2) -- check Library
#     Column  12 =  c_s_area  = Area exposed by the cluster excluding the interface (in Angstrom^2)
#     Column  13 =   Esurf    = Exposed surface energy (in J/m^2) of the cluster (not interface) -- check Library
#     Column  14 =    Ecoh    = Cohesion energy per cluster atom (in eV/atom) == (Ecluster -( N * Eatom))/N
#     Column  15 =    Eadh    = Cluster adhesion energy (in eV) == Esystem - (Esurface + Ecluster)
#     Column  16 =     Eb     = Binding energy per cluster atom (in eV/atom) == (Esystem -(Esurface + N * Eatom))/N
#     Column  17 =   Etotal   = Total energy of the system (in eV)
#     Column  18 = structure_path

	un_labels.pop(0)
	site_a = un_labels[2][3:]
	site_b = un_labels[3][3:]
	labels = ["$NP_{n}$", "# $NP^{interface}$", "# " + site_a + "$^{interface}$", "# " + site_b + "$^{interface}$",
			  "$\overline{coordination_{NP^{interface}}}$",
			  "$\overline{d_{NP-" + site_a + "}^{min}}$ $(\AA)$", "$\overline{d_{NP-" + site_b + "}^{min}}$ $(\AA)$",
			  "$\overline{d^{interface}}$ $(\AA)$", "$d_{ \overline{support} - NP_{CoM}}$ $(\AA)$",
			  "$\overline{GCN_{NP}}$", "$Area_{NP^{interface}}$ $(\AA^{2})$", "$Area_{NP^{external}}$ $(\AA^{2})$",
			  "$\gamma_{NP}$ $(J \cdot m^{\minus 2})$", "$E_{Coh}$ $(eV \cdot atom^{\minus 1})$",
			  "$E_{Adh}$ $(eV)$", "$E_{B}$ $(eV \cdot atom^{\minus 1})$", "$E_{Total}$ $(eV)$"]

	data = []
	system_name = []
	symbol = []
	for i in range(len(dataset)):
		energies = [float(dataset[i][14]), float(dataset[i][15]), float(dataset[i][16]), float(dataset[i][17])]
		mask = [(dataset[i][-1], labels[n+13], e) for n, e in enumerate(energies) if e > 0.25]
		if len(mask) > 0:
			[print(m) for m in mask]
		else:
			data.append([float(j)for j in dataset[i][:-2]])
			system_name.append(str(dataset[i][-1]))
			symbol.append(str(dataset[i][-1].split("/")[1] + dataset[i][-1].split("/")[2]))

	return data, labels, symbol, labels_line, system_name


def Display(xlabel, ylabel, xlim, ylim, trend_label):
	ax.set_xlabel(xlabel, fontsize=14)
	ax.set_ylabel(ylabel, fontsize=14)
	ax.tick_params(axis='both', labelsize=12)
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	if trend_label != "":
		ax.set_title(trend_label, position=(0.5, 0.95))
# make BLUE the structures used for fitting, e.g. 2, 3, 4a, 5a, ...
	fitted_structures = ["2", "3", "4a", "5a", "6a", "7", "8", "9", "10", "11"]
	legend = ax.legend(loc='center left', bbox_to_anchor=(1.04, 0.5), borderaxespad=0.) #loc='best')
	for text in legend.get_texts():
		text_start = text.get_text().split("$")[0][2:]
		if text_start in fitted_structures:
			plt.setp(text, color="b")
		else:
			plt.setp(text, color="k")
	figure.tight_layout()#rect=[0, 0, 0.75, 1])
	plt.ion()
	plt.show()
	SaveFig()
	figure.clf()


def SaveFig():
	answer = str(input("Would you like to save the figure (y/n)?\n"))
	if answer == "y":
		figure_out_name = str(input("What would it be the figure name (a word & no format)?\n"))
		plt.savefig(figure_out_name + ".svg", figsize=(14, 16), clear=True,
										bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)


def Validation(name, x, y, imarker, icolour):
	deviation = [np.abs((y[i] - x[i])/x[i]) for i in range(len(x))]
# Add label to each point if deviation is larger than 0.5
	if max(deviation) > 1.0:
		for i in range(len(x)):
			if x[i] <= 0 and np.abs(y[i] - x[i]) > 0.3:
				ax.text(x[i]+0.02, y[i]+0.02, str(i+1))
	ax.plot(x, y,  marker=imarker, color=icolour, linestyle="None", alpha=0.5,
			 label=str(name) + "$\cdot \\tau \leq$ " + str(round(max(deviation), 1)) + "%")
	return max(deviation)

########################################################################################################################
data_a, labels_a, symbol_a, labels_line_a, system_name_a = get_data(sys.argv[1])
data_b, labels_b, symbol_b, labels_line_b, system_name_b = get_data(sys.argv[2])

# CONTROL of comparable files
if len(labels_a) != len(labels_b):
	print(" * The number of columns is {:d} and {:d} for {} and {} *".format(len(labels_a), len(labels_b), *sys.argv[1:]))
	exit()
if len(data_a) != len(data_b):
	print(" ** The number of rows is {:d} and {:d} for {} and {} **".format(len(data_a), len(data_b), *sys.argv[1:]))
	data_missing = []
	for i in range(len(system_name_a)):
		if system_name_a[i] not in system_name_b:
			data_missing.append([i + labels_line_a + 3, system_name_a[i+1], sys.argv[2]])
	for i in range(len(system_name_b)):
		if system_name_b[i] not in system_name_a:
			data_missing.append([i + labels_line_b + 3, system_name_b[i+1], sys.argv[1]])
	if len(data_missing) == 0:
		print("\tThere are no systems missing")
	else:
		for i in data_missing:
			print("\tIn line {:d}, the system {} is missing in {}".format(*i))
	exit()
for i in range(len(symbol_a)):
	if symbol_a[i] != symbol_b[i]:
		print(" *** The entries in row {:d} are different in for {} and {} ***".format(i, *sys.argv[1:]))
		exit()
print("   The dataset size is: ", len(system_name_a))

#--------------------------------------- Validation ---------------------------------------
trend_file = open("Validation_Summary.txt", 'w+')
names = []
for i in symbol_a:
	if i not in names:
		names.append(i)
for n in range(13, len(labels_a)):
	figure = plt.figure(figsize=(8, 6), clear=True)
	ax = figure.add_subplot()
	axis = []
	for i, name in enumerate(names):
		n_m = i
		n_c = i
		if i >= 2*len(icolour):
			n_c = i - 2*len(icolour)
		elif i >= len(icolour):
			n_c = i - len(icolour)
		if n >= len(imarker):
			n_m = i - len(imarker)

		x = [data_a[j][n+1] for j in range(len(data_a)) if symbol_a[j] == name]
		y = [data_b[j][n+1] for j in range(len(data_b)) if symbol_b[j] == name]
		axis += x + y
		max_deviation = Validation(name, x, y, imarker[n_m], icolour[n_c])
		trend_file.write("# Column {}: {}\tSystem {}\tMaximum Absolute Error: \u03C4\u2264{:<1.2f}\n"
						 .format(n, labels_a[n], name, max_deviation))

	axis_max = max(axis) + np.abs(max(axis)*0.05)
	axis_min = min(axis) - np.abs(min(axis)*0.05)
	if axis_max - axis_min <= 1:
		axis_max = axis_max + np.abs(axis_min*0.5)
		axis_min = axis_min - np.abs(axis_min*0.5)
#	ax.text((axis_max - axis_min)*1.8/3, (axis_max - axis_min)*0.05, "Samples = " + str(len(symbol_a)))
	ax.plot([axis_min, axis_max], [axis_min, axis_max], "k-", lw=1.5)
	Display(labels_a[n], "Predicted " + labels_a[n], [axis_min, axis_max], [axis_min, axis_max], "Samples = " + str(len(symbol_a)))
trend_file.close()
