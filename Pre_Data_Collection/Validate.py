'''

	USAGE: ~.py input.dat
   	input: average_coordination E_Coh (eV.atom^-1) >> # Symbols Path

'''

import sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

icolour = ["b", "r", "k", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] ## n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))]


def get_data(data):
	ic = []									# contains the number of interface_cluster_atoms
	icc = []								# contains the average coordination of the cluster atoms at the interface
	id = []									# contains the average distance from the cluster interface atoms to the surface neighbouring atoms
	isd_a = []								# contains the average of the shortest distances from the  interfacial atoms in the cluster to the preferable surface site[0]
	isd_b = []								# contains the average of the shortest distances from the  interfacial atoms in the cluster to the preferable surface site[1]
	e = []								# contains the DFT calculated adhesion energies
	predicted_e = []
	name = []
	for i in range(len(data)):
		if float(data[i][6]) <= 0.25:
			ic.append(float(data[i][0])+float(data[i][1]))
			icc.append(float(data[i][2]))
			id.append(float(data[i][3]))
			isd_a.append(float(data[i][4]))
			isd_b.append(float(data[i][5]))
			e.append(float(data[i][6]))
			predicted_e.append(float(data[i][7]))
			name.append(data[i][-1].split("/")[-1])

	return ic, icc, id, isd_a, isd_b, e, predicted_e, name


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
		plt.savefig(figure_out_name + ".svg", figsize=(12, 10), clear=True,
										bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)


def trend_lineal(x, y, symbol, xlim, colour, marker, line):
	popt = np.polyfit(x, y, 1)
	p = np.poly1d(popt)
	r2 = 1-np.sqrt(sum([(y[i] - p(y[i]))**2 for i in range(len(y))])/sum(i*i for i in y))
	a, b = popt
	trend_label = "Predicted_EAdh= {:.5f} * E_Adh + {:.5f}".format(a, b)
	plt.plot(np.linspace(xlim[0], xlim[1], 150), p(np.linspace(xlim[0], xlim[1], 150)), color=colour, linestyle=line)
#	plt.plot(x, y, marker=marker, color=colour, linestyle="None", label=str(symbol) + "$\cdot R^{2}$= "+str(round(r2, 2)))

	return trend_label, r2


def Validation_3D(ele, x, y, name, imarker, icolour):

	trend_label, r2 = trend_lineal(x, y, ele, [min(x), max(x)], "g", "*", "-")

	deviation = [(y[i] - x[i]) for i in range(len(x))]
# Add label to each point
	trend_file = open("Deviated_Addresses_EAdh.txt", 'w+')
	trend_file.write("# Data that falls in a margin of error > 0.2\n")
	for i in range(len(x)):
		if x[i] <= 0 and np.abs(deviation[i]) > 0.2:
			plt.text(x[i]+0.02, y[i]+0.02, str(i+1))
			trend_file.write("{} --- {}\n".format(i+1, name[i]))
	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None", label=str(ele) + "$\cdot \\tau \leq$ " +\
																			   str(np.abs(round(max(deviation), 1))))
	trend_file.close()
	return max(deviation), trend_label, r2

########################################################################################################################
dataset_length = 0
symbol = []
name = {}
ic = {}
icc = {}
id = {}
isd_a = {}
isd_b = {}
e = {}
predicted_e = {}
reference_e = {}
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	dataset_length += len(data)
	print("   The dataset size is: ", len(data))
	symbol.append(data[0][-1].split("/")[2]) # [0])					# contains the list of systems' name
	ic[symbol[-1]], icc[symbol[-1]], id[symbol[-1]], isd_a[symbol[-1]], isd_b[symbol[-1]], e[symbol[-1]],\
	predicted_e[symbol[-1]], name[symbol[-1]] = get_data(data)
z_min = min([min(e[sym]) for sym in symbol]) - np.abs(min([min(e[sym]) for sym in symbol]))*0.1
print("   The TOTAL dataset size is: ", dataset_length)

#--------------------------------------- Validation ---------------------------------------
trend_file = open("Predicted_EAdh.txt", 'w+')
for n, sym in enumerate(symbol):
	max_deviation, trend_label, r2 = Validation_3D(sym, e[sym], predicted_e[sym], name[sym], imarker[n], icolour[n])
	trend_file.write("# System = {}\n#\tE_Adh (eV)\t\tMaximum Absolute Error: \u03C4\u2264{:<1.2f} eV\n".format(sym,
																								np.abs(max_deviation)))
	trend_file.write("\t{}\t\tR\u00b2={:<1.2f}\n".format(trend_label, round(r2, 2)))

plt.plot([z_min, 0], [z_min, 0], "k-", lw=1.5)
Display("$E_{Adh}$ (eV)", "Predicted $E_{Adh}$ (eV)", [z_min, 0], [z_min, 0], "")
trend_file.close()
