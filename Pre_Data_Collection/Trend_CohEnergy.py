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


# NOTE that 12 is the bulk coordination for FCC metals
# Coh is NORMALISED by the bulk Coh
def logarithm(x, a):#, b, c): #, d):
	return np.log(a)/np.log(a/(a+bulk_coord)) - (1/np.log(a/(a+bulk_coord)))*np.log(a+x)
#	return 1/np.log(a/(a+b)) * (np.log(a) - c* np.log(a+x))
#	d/np.log(a/(a+b)) * (np.log(a) - c * np.log(a+x))


def cubic(x, a, b, c, d):
	return a*x**3 + b*x**2 + c*x + d


def trend_function(x, y, bulk_coord, symbol, colour, marker, line, name):
	plt.plot([0, bulk_coord], [0, 1], marker=marker, color=colour, fillstyle="none", linestyle="None", markersize=3)
	weights = list(np.ones(len(x)))
	x.append(0.)
	y.append(0.)
	name.append("")
	weights.append(0.1)
	x.append(bulk_coord)
	y.append(1)
	name.append("bulk")
	weights.append(0.1)
#	popt, pcov = curve_fit(cubic, x, y, sigma=weights)
	popt, pcov = curve_fit(logarithm, x, y, sigma=weights)
#	r2 = 1-np.sqrt(sum([(y[i] - cubic(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
	r2 = 1-np.sqrt(sum([(y[i] - logarithm(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
#	a, b, c, d = popt
#	trend_label = "Cubic: ax^3 + bx^2 +  cx + d = {:.2f}, {:.2f}, {:.2f}, {:.2f},".format(a, b, c, d)
	trend_label = "logarithmic: {:.2f}".format(*popt)
	print(trend_label)
#	plt.plot(np.linspace(0, 12, 150), cubic(np.linspace(0, 12, 150), *popt), color=colour, linestyle=line)
	plt.plot(np.linspace(0, 12, 150), logarithm(np.linspace(0, 12, 150), *popt), color=colour, linestyle=line)

	magic_clusters = ["116", "13", "140", "146", "147", "19", "201", "225", "38", "44", "55", "6a", "79", "85"]
	magic_x = [x[i] for i in range(len(x)) if str(name[i][2:]) in magic_clusters]
	magic_y = [y[i] for i in range(len(y)) if str(name[i][2:]) in magic_clusters]
	new_x = [x[i] for i in range(len(x)) if str(name[i][2:]) not in magic_clusters]
	new_y = [y[i] for i in range(len(y)) if str(name[i][2:]) not in magic_clusters]
	plt.plot(magic_x, magic_y, marker=marker, color=colour, fillstyle='none', linestyle="None", markersize=3)
	plt.plot(new_x, new_y, marker=marker, color=colour, linestyle="None", markersize=3,
			 label=str(symbol) + "$\cdot R^{2}$= "+str(round(r2, 2)))
# Add names next to the point
#	for i in range(len(name[:-2])):
#		plt.text(x[i]+0.05, y[i]+0.01, str(name[i]), color=colour)

	return trend_label, popt, r2


def Validation(ele, coord, coh_e, popt, imarker, icolour):
	x = coh_e
	y = [logarithm(i, *popt) for i in coord]          		# in eV.atom^-1
	max_deviation = max([np.abs(y[i] - x[i]) for i in range(len(x))])
	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None", markersize=3,
			 label=str(ele) + "$\cdot \\tau \leq$ " + str(round(max_deviation, 1)) + "eV")

	return max_deviation


########################################################################################################################
trend = {}
r = {}
name = []
y_min = 1
y_max = 1
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol = data[0][-1].split("/")[0]					# contains the list of systems' name
	names = [data[i][-2] for i in range(len(data))]
	coord, coh_e, coh_e_bulk, bulk_coord = get_data(data)
	trend_label, trend[symbol], r[symbol] = trend_function(coord, coh_e, bulk_coord, symbol, icolour[n-1], imarker[n-1], iliner[n], names)
Display("$\overline{coordination_{NP}}$", "$\\frac{E_{Coh}}{E_{Coh}^{Bulk}}$", [-0.15, 12.15], [-0.02, 1.02], trend_label)
########################################################################################################################
trend_file = open("Interpolation_ECoh_Trends.txt", 'w+')
plt.plot([-0.02, 1.02], [-0.02, 1.02], "k-", lw=1.5)
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol = data[0][-1].split("/")[0]					# contains the list of systems' name
	coord, coh_e, coh_e_bulk, bulk_coord = get_data(data)
	max_deviation = Validation(symbol[:2], coord, coh_e, trend[symbol], imarker[n-1], icolour[n-1])
	trend_file.write("# E_Coh/E_Coh^Bulk\tLogarithm interpolation: 1/E_Coh^Bulk * log(a)/np.log(a/(a+bulk_coord)) - (1/np.log(a/(a+bulk_coord)))*np.log(a+cc)\n")
	trend_file.write("{}\ta={:<3.5f}\tR\u00B2={:<1.2f}\t\u03C4\u2264{:<1.1f}eV\n"
					 .format(symbol, trend[symbol][0], r[symbol], np.abs(round(max_deviation, 1))))
#	trend_file.write("{}\ta={:<3.5f}\tb={:<5.5f}\tc={:<3.5f}\tR\u00B2={:<1.2f}\t\u03C4\u2264{:<1.1f}eV\n"
#					 .format(symbol, trend[symbol][0], trend[symbol][1], trend[symbol][2], r[symbol], np.abs(round(max_deviation, 1))))
Display("$\\frac{E_{Coh}}{E_{Coh}^{Bulk}}$", "Predicted $\\frac{E_{Coh}}{E_{Coh}^{Bulk}}$", [-0.02, 1.02], [-0.02, 1.02], "")
