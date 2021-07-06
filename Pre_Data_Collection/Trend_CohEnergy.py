'''

	USAGE: ~.py input.dat
   	input: average_coordination E_Coh (eV.atom^-1) >> # Symbols Path

'''

import sys
import numpy as np
from scipy.optimize import minimize, curve_fit
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt



icolour = ["b", "r", "k", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] ## n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))]


def get_data(data):
	coord = [float(data[i][0]) for i in range(len(data))]			# contains the average coordinations
	coh_e_bulk = sum([float(data[i][3]) for i in range(len(data))])/(len(data))	# contains the average bulk cohesion energies in eV.atom^-1
	coh_e = [float(data[i][1])/coh_e_bulk for i in range(len(data))]			# contains the (cohesion energies) / (Ecoh bulk average)
	bulk_coord = sum([float(data[i][2]) for i in range(len(data))])/len(data)	# average coord in the bulk (fcc~12)

	return coord, coh_e, coh_e_bulk, bulk_coord


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


def logarithm(x, a, b, c):
	return np.log(a)/np.log(a/(a+b)) - (c/np.log(a/(a+b)))*np.log(a+x)


def trend_logarithm(x, y, bulk_coord, symbol, colour, marker, line):
	plt.plot(x, y, marker=marker, color=colour, linestyle="None")
	plt.plot([0, bulk_coord], [0, 1], marker=marker, color=colour, fillstyle="none", linestyle="None")
	weights = list(np.ones(len(x)))

	x.append(0.)
	x.append(bulk_coord)
	weights.append(0.5)
	y.append(0.)
	y.append(1)
	weights.append(0.5)

	popt, pcov = curve_fit(logarithm, x, y, sigma=weights)
	r2 = 1-np.sqrt(sum([(y[i] - logarithm(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
	a, b, c = popt
	trend_label = "Logarithm: a, b, c =" + str(round(a, 5)) + ", " + str(round(b, 5)) + ", " + str(round(c, 5))
	plt.plot(np.linspace(0, 12, 150), logarithm(np.linspace(0, 12, 150), *popt), color=colour, linestyle=line,
			 label=str(symbol) + "$\cdot R^{2}$= "+str(round(r2, 2)))

	return trend_label, popt, r2


def Validation(ele, coord, coh_e, popt, imarker, icolour):
	x = coh_e
	y = [trend_logarithm(i, *popt) for i in coord]          		# in eV.atom^-1
	max_deviation = max([(y[i] - x[i])*100/x[i] for i in range(len(x))])

	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None",
			 label=str(ele) + "$\cdot \\tau \leq$ " + str(np.abs(round(max_deviation, 1))) + "%")

	return max_deviation


########################################################################################################################
trend = {}
r = {}
y_min = 1
y_max = 1
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol = data[0][-2].split("/")[0]					# contains the list of systems' name

	coord, coh_e, coh_e_bulk, bulk_coord = get_data(data)
	if y_min > min(coh_e):
		y_min = min(coh_e)
	if y_max < max(coh_e):
		y_max = max(coh_e)

	trend_label, trend[symbol], r[symbol] = trend_logarithm(coord, coh_e, bulk_coord, symbol[:2], icolour[n-1], imarker[n-1], iliner[n])
Display("Average Coordination", "$\\frac{E_{Coh}}{E_{Coh}^{Bulk}}$", [-0.15, 12.15], [-0.02, 1.02], trend_label)

trend_file = open("Interpolation_ECoh_Trends.txt", 'w+')
plt.plot([y_min, y_max], [y_min, y_max], "k-", lw=1.5)
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol = data[0][-2].split("/")[0]					# contains the list of systems' name

np.log(a)/np.log(a/(a+b)) - (c/np.log(a/(a+b)))*np.log(a+x)

	coord, coh_e, coh_e_bulk, bulk_coord = get_data(data)
	max_deviation = Validation(symbol[:2], coord, coh_e, trend[symbol], imarker[n], icolour[n], y_min, y_max)
		trend_file.write("# E_Coh/E_Coh^Bulk\tLogarithm interpolation: 1/E_Coh^Bulk * a - (b + x**c)/d\n")
		trend_file.write("\ta={:<3.5f}\tb={:<5.5f}\tc={:<3.5f}\td={:<5.5f}\tR\u00B2={:<1.2f}\t\u03C4\u2264{:<1.1f}%\n" .format(
			area_popt[ele][0], area_popt[ele][1], area_popt[ele][2], area_popt[ele][3], r2, np.abs(round(max_deviation, 1))))

Display("Area ($\\AA ^{2}$)", "Predicted Area ($\\AA ^{2}$)", [x_min, x_max], [x_min, x_max], set(elements))
