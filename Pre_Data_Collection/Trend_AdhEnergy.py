'''

	USAGE: ~.py input.dat
   	input: average_coordination E_Coh (eV.atom^-1) >> # Symbols Path

'''

import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt



icolour = ["b", "r", "k", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] ## n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))]


def get_data(data):
	ic = [float(data[i][0]) for i in range(len(data))]			# contains the number of interface_cluster_atoms
	icc = [float(data[i][1]) for i in range(len(data))]			# contains the average coordination of the cluster atoms at the interface
	id = [float(data[i][2]) for i in range(len(data))]			# contains the average distance from the cluster interface atoms to the surface neighbouring atoms
	idr = [float(data[i][3]) for i in range(len(data))]			# contains the shortest distance from the reference (FIRST interfacial atom in the cluster) to the a surface neighbour
	adh_e = [float(data[i][4]) for i in range(len(data))]		# contains the DFT calculated adhesion energies

	return ic, icc, id, idr, adh_e


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


def logarithm(x, a, b, c, d):
	return d/np.log(a/(a+b)) * (np.log(a) - c * np.log(a+x))


def cubic(x, a, b, c, d):
	return a*x**3 + b*x**2 + c*x + d


def trend_function(x, y, symbol, colour, marker, line):
	popt, pcov = curve_fit(cubic, x, y)
	r2 = 1-np.sqrt(sum([(y[i] - cubic(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
	a, b, c, d = popt
	trend_label = "Cubic: ax^3 + bx^2 +  cx + d = {:.2f}, {:.2f}, {:.2f}, {:.2f},".format(a, b, c, d)
	plt.plot(np.linspace(0, 12, 150), cubic(np.linspace(0, 12, 150), *popt), color=colour, linestyle=line)
	plt.plot(x, y, marker=marker, color=colour, linestyle="None", label=str(symbol) + "$\cdot R^{2}$= "+str(round(r2, 2)))

	return trend_label, popt, r2


def Validation(ele, coord, coh_e, popt, imarker, icolour):
	x = coh_e
	y = [cubic(i, *popt) for i in coord]          		# in eV.atom^-1
	max_deviation = max([(y[i] - x[i])*100/x[i] for i in range(len(x))])
	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None",
			 label=str(ele) + "$\cdot \\tau \leq$ " + str(np.abs(round(max_deviation, 1))) + "%")

	return max_deviation


########################################################################################################################
trend = {}
r = {}
y_min = 1
y_max = 1
x_min = 1
x_max = 1
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol = data[0][-1].split("/")[0]					# contains the list of systems' name
	ic, icc, id, idr, adh_e = get_data(data)
	if x_min > min(id):
		x_min = min(id)
	if x_max < max(id):
		x_max = max(id)
	if y_min > min(adh_e):
		y_min = min(adh_e)
	if y_max < max(adh_e):
		y_max = max(adh_e)
#------------------------------------- Interfacial Distance ------------------------
	trend_label, trend[symbol], r[symbol] = trend_function(id, adh_e, symbol, icolour[n-1], imarker[n-1], iliner[n])
x_min = x_min*0.8
x_max = x_max*1.5
y_min = y_min - np.abs(y_min*0.1)
y_max = y_max*1.05
Display("Interface Average Distance ($\\AA$)", "$E_{Adh}$ (eV)", [x_min, x_max], [y_min, y_max], trend_label)

trend_file = open("Interpolation_EAdh_Trends.txt", 'w+')
y_min = y_min*0.8
y_max = y_max*1.05
plt.plot([y_min, y_max], [y_min, y_max], "k-", lw=1.5)
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol = data[0][-1].split("/")[0]					# contains the list of systems' name
	ic, icc, id, idr, adh_e = get_data(data)
	max_deviation = Validation(symbol[:2], id, adh_e, trend[symbol], imarker[n-1], icolour[n-1])
	trend_file.write("# E_Coh/E_Coh^Bulk\tLogarithm interpolation: 1/E_Coh^Bulk * 1/log(a/(a+b)) * (log(a) - c* log(a+cc))\n")
	trend_file.write("{}\ta={:<3.5f}\tb={:<5.5f}\tc={:<3.5f}\tR\u00B2={:<1.2f}\t\u03C4\u2264{:<1.1f}%\n"
					 .format(symbol, trend[symbol][0], trend[symbol][1], trend[symbol][2], r[symbol], np.abs(round(max_deviation, 1))))
Display("$\\frac{E_{Coh}}{E_{Coh}^{Bulk}}$", "Predicted $\\frac{E_{Coh}}{E_{Coh}^{Bulk}}$", [y_min, y_max], [y_min, y_max], "")
