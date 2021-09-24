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
	labels_line = len([1 for i in range(len(ifile)) if ifile[i].startswith("#") is True]) - 2
#     Column   0 =     N      = Total number of atoms forming the cluster
#     Column   1 =    i_c     = Number of cluster atoms coordinating the surface
#     Column   2 =    cs_O    = Number of surface sites coordinating with the cluster
#     Column   3 =   cs_Mg    = Number of surface sites coordinating with the cluster
#     Column   4 =    i_cc    = Average coordination of cluster atoms at the interface within the cluster only
#     Column   5 =     cc     = Average atomic coordination within the cluster
#     Column   6 =   dist_O   = Average of minimum distances (in Å) between the surface sites and the clusters atoms
#     Column   7 =  dist_Mg   = Average of minimum distances (in Å) between the surface sites and the clusters atoms
#     Column   8 = cs_height  = Number of surface sites coordinating with the cluster
#     Column   9 =   Zdist    = Distance (in Å) between the average surface hight and the cluster's centre of mass
#     Column  10 =    GCN     = Average generalised coordination number for the atoms in the cluster excluding the coordination with the support
#     Column  11 =  c_i_area  = Cluster interface area (in Angstrom^2) -- check Library
#     Column  12 =  c_s_area  = Area exposed by the cluster excluding the interface (in Angstrom^2)
#     Column  13 =   Esurf    = Exposed surface energy (in J/m^2) of the cluster (not interface) -- check Library
#     Column  14 =    Ecoh    = Cohesion energy per cluster atom (in eV/atom) == (Ecluster -( N * Eatom))/N
#     Column  15 =    Eadh    = Cluster adhesion energy (in eV) == Esystem - (Esurface + Ecluster)
#     Column  16 =     Eb     = Binding energy per cluster atom (in eV/atom) == (Esystem -(Esurface + N * Eatom))/N
#     Column  17 =   Etotal   = Total energy of the system (in eV)
#     Column  18 = structure_path

	labels = []
	for i in range(len(ifile[labels_line].split())):
		if i in [6, 7, 9]:
			labels.append(ifile[labels_line].split()[i] + " $(\AA)$")
		elif i in [11, 12]:
			labels.append(ifile[labels_line].split()[i] + " $(\AA^{2})$")
		elif i in [13]:
			labels.append(ifile[labels_line].split()[i] + " $(J \cdot m^{\minus 2})$")
		elif i in [14, 16]:
			labels.append(ifile[labels_line].split()[i] + " $(eV \cdot atom^{\minus 1})$")
		elif i in [15, 17]:
			labels.append(ifile[labels_line].split()[i] + " $(eV)$")
		else:
			labels.append(ifile[labels_line].split()[i])

	data = []
	symbol = []
	for i in range(len(dataset)):
		if float(dataset[i][14]) and float(dataset[i][15]) and float(dataset[i][16]) and float(dataset[i][17]) <= 0.25 :
			data.append([float(j) for j in dataset[i][:-2]])
			symbol.append(str(dataset[i][-1].split("/")[1] + dataset[i][-1].split("/")[2]))

	return data, labels, symbol


def Display(xlabel, ylabel, xlim, ylim, trend_label):
	plt.xlabel(str(xlabel), fontsize=14)
	plt.ylabel(str(ylabel), fontsize=14)
	plt.tick_params(axis='both', labelrotation=0, labelsize=12)               # custimise tick labels
#	plt.xticks(np.arange(int(xlim[0]), int(xlim[1]), 1))   # Xmin,Xmax,Xstep
	plt.xlim(xlim)
	plt.ylim(ylim)
	if trend_label != "":
		plt.title(trend_label)
	plt.legend(loc='best')
	plt.tight_layout()
	plt.ion()
	plt.show()
	SaveFig()
	plt.clf()


def SaveFig():
	answer = str(input("Would you like to save the figure (y/n)?\n"))
	if answer == "y":
		figure_out_name = str(input("What would it be the figure name (a word & no format)?\n"))
		plt.savefig(figure_out_name + ".svg", figsize=(16, 16), clear=True,
										bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)


def Validation(name, x, y, imarker, icolour):
	deviation = [np.abs(y[i] - x[i]) for i in range(len(x))]
# Add label to each point if deviation is larger than 0.5
	if max(deviation) > 0.5:
		for i in range(len(x)):
			if x[i] <= 0 and np.abs(y[i] - x[i]) > 0.3:
				plt.text(x[i]+0.02, y[i]+0.02, str(i+1))
	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None", alpha=0.5,
			 label=str(name) + "$\cdot \\tau \leq$ " + str(round(max(deviation), 1)))
	return max(deviation)

########################################################################################################################
data_a, labels_a, symbol_a = get_data(sys.argv[1])
data_b, labels_b, symbol_b = get_data(sys.argv[2])

# CONTROL of comparable files
if len(labels_a) != len(labels_b):
	print(" * The number of columns is {:d} and {:d} for {} and {} *".format(len(labels_a), len(labels_b), *sys.argv[1:]))
	exit()
if len(data_a) != len(data_b):
	print(" ** The number of rows is {:d} and {:d} for {} and {} **".format(len(data_a), len(data_b), *sys.argv[1:]))
	exit()
for i in range(len(data_a)):
	if symbol_a[i] != symbol_b[i]:
		print(" *** The entries in row {:d} are different in for {} and {} ***".format(i, *sys.argv[1:]))
		exit()
print("   The dataset size is: ", len(symbol_a))

#--------------------------------------- Validation ---------------------------------------
trend_file = open("Validation_Summary.txt", 'w+')
name = symbol_a[0]
for n, label in enumerate(labels_a):
	i = 0
	if symbol_a[i] != name:
		name = symbol_a[i]
		i += 1
	n_m = i
	n_c = i
	if i >= 2*len(icolour):
		n_c = i - 2*len(icolour)
	elif i >= len(icolour):
		n_c = i - len(icolour)
	if n >= len(imarker):
		n_m = i - len(imarker)

	max_deviation = Validation(name, data_a[:][n], data_b[:][n], imarker[n_m], icolour[n_c])
	trend_file.write("# Column {}: {}\tSystem {}\tMaximum Absolute Error: \u03C4\u2264{:<1.2f}\n"
					 .format(n, label, name, max_deviation))
	axis_max = max(data_a[:][n] + data_b[:][n])
	axis_min = min(data_a[:][n] + data_b[:][n])
	plt.text((axis_max - axis_min)*1.8/3, (axis_max - axis_min)*0.05, "Samples = " + str(len(symbol_a)))
	plt.plot([axis_min, axis_max], [axis_min, axis_max], "k-", lw=1.5)
	Display(label, "Predicted " + label, [axis_min, axis_max], [axis_min, axis_max], "")
trend_file.close()
