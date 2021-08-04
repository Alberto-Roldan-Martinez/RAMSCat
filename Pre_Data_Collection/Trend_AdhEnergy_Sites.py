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
	ic = []									# contains the number of interface_cluster_atoms
	icc = []								# contains the average coordination of the cluster atoms at the interface
	id = []									# contains the average distance from the cluster interface atoms to the surface neighbouring atoms
	isd_a = []								# contains the average of the shortest distances from the  interfacial atoms in the cluster to the preferable surface site[0]
	isd_b = []								# contains the average of the shortest distances from the  interfacial atoms in the cluster to the preferable surface site[1]
	adh_e = []								# contains the DFT calculated adhesion energies
	scaled_adh_e = []
	for i in range(len(data)):
		if float(data[i][0]) * float(data[i][1]) > 0:
			if float(data[i][5]) < 0:
				ic.append(float(data[i][0]))
				icc.append(float(data[i][1]))
				id.append(float(data[i][2]))
				isd_a.append(float(data[i][3]))
				isd_b.append(float(data[i][4]))
				adh_e.append(float(data[i][5]))
				scaled_adh_e.append(float(data[i][5])/float(data[i][0])) # (float(data[i][0]) * float(data[i][1])))

	return ic, icc, id, isd_a, isd_b, adh_e, scaled_adh_e


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

def morse(x, a, d_eq, r_eq):
	return d_eq * (np.exp(-2*a*(x - r_eq)) - 2 * np.exp(-a*(x - r_eq)))     # MORSE potential

def trend_morse(x, y, symbol, xlim, colour, marker, line):
	popt, pcov = curve_fit(morse, x, y, bounds=([0., 0.001, 0.1], [10, 10, 10]))
	r2 = 1-np.sqrt(sum([(y[i] - morse(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
	a, b, c = popt
	trend_label = "Morse: {:.2f} * (exp(-2*{:.2f}*(id - {:.2f})) - 2 * exp(-{:.2f}*(id - {:.2f})))".format(b, a, c, a, c)
	plt.plot(np.linspace(xlim[0], xlim[1], 150), morse(np.linspace(xlim[0], xlim[1], 150), *popt), color=colour, linestyle=line)
	plt.plot(x, y, marker=marker, color=colour, linestyle="None", label=str(symbol) + "$\cdot R^{2}$= "+str(round(r2, 2)))

	return trend_label, popt, r2


def Validation(ele, x0, z0, popt, imarker, icolour):
	x = z0
	y = morse(x0, *popt)
	max_deviation = max([(y[i] - x[i])*100/x[i] for i in range(len(x))])
	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None",
			 label=str(ele) + "$\cdot \\tau \leq$ " + str(np.abs(round(max_deviation, 1))) + "%")

	return max_deviation


########################################################################################################################
symbol = []
ic = {}
icc = {}
id = {}
isd_a = {}
isd_b = {}
adh_e = {}
scaled_adh_e = {}
id_trend = {}
a_trend = {}
b_trend = {}
id_r = {}
a_r = {}
b_r = {}
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol.append(data[0][-1].split("/")[2]) # [0])					# contains the list of systems' name
	ic[symbol[-1]], icc[symbol[-1]], id[symbol[-1]], isd_a[symbol[-1]], isd_b[symbol[-1]], adh_e[symbol[-1]], scaled_adh_e[symbol[-1]] = get_data(data)
dx_min = min([min(id[sym]) for sym in symbol])*0.9
dx_max = max([max(id[sym]) for sym in symbol])*1.3
x_min = min([min(isd_a[sym]) for sym in symbol])*0.9
x_max = max([max(isd_a[sym]) for sym in symbol])*1.3
y_min = min([min(isd_b[sym]) for sym in symbol])*0.9
y_max = max([max(isd_b[sym]) for sym in symbol])*1.3
z_min = min([min(scaled_adh_e[sym]) for sym in symbol]) - np.abs(min([min(scaled_adh_e[sym]) for sym in symbol]))*0.1
z_max = max([max(scaled_adh_e[sym]) for sym in symbol])*1.3
#-------------------------------------  Interface Distance  ------------------------
for n, sym in enumerate(symbol):
	trend_label, id_trend[sym], id_r[sym] = trend_morse(id[sym], scaled_adh_e[sym], sym, [dx_min, dx_max], icolour[n], imarker[n], iliner[n+1])
Display("id ($\\AA$)", "$E_{Adh}^{Scaled}$ $(eV \cdot atom^{-1})$", [dx_min, dx_max], [z_min, 0], trend_label)
#--------------------------------------- Validation ---------------------------------------
trend_file = open("Interpolation_EAdh_Trends.txt", 'w+')
for n, sym in enumerate(symbol):
	max_deviation = Validation(sym, id[sym], scaled_adh_e[sym], id_trend[sym], imarker[n], icolour[n])
	a, d_eq, r_eq = id_trend[sym]
	trend_file.write("# Scaled_E_Adh (eV . n_interface_cluster_atoms^-1)\n#  Interface Distance\tMorse interpolation:"
					 " d_eq * (exp(-2 * a * (D - r_eq)) - 2 * exp(-a * (D - r_eq)))\n")
	trend_file.write("{}\ta={:<5.5f}\td_eq={:<5.5f}\tr_eq={:<5.5f}\tR\u00B2={:<1.2f}\t\u03C4\u2264{:<1.1f}%\n"
					 .format(sym, a, d_eq, r_eq, id_r[sym], np.abs(max_deviation)))
plt.plot([z_min, 0], [z_min, 0], "k-", lw=1.5)
Display("$E_{Adh}$ (eV)", "Predicted $E_{Adh}$ (eV)", [z_min, 0], [z_min, 0], "")
trend_file.write("\n")
trend_file.close()

#-------------------------------------  Distance A ------------------------
for n, sym in enumerate(symbol):
	trend_label, a_trend[sym], a_r[sym] = trend_morse(isd_a[sym], scaled_adh_e[sym], sym, [x_min, x_max], icolour[n], imarker[n], iliner[n+1])
#x_min = min([a_trend[sym][2] for sym in symbol])*0.85
#z_min = -1*min([a_trend[sym][1] for sym in symbol]) - np.abs(min([a_trend[sym][1] for sym in symbol]))*0.1

Display("isd_a ($\\AA$)", "$E_{Adh}^{Scaled}$ $(eV \cdot atom^{-1})$", [x_min, x_max], [z_min, 0], trend_label)
#--------------------------------------- Validation ---------------------------------------
trend_file = open("Interpolation_EAdh_Trends.txt", 'a+')
for n, sym in enumerate(symbol):
	max_deviation = Validation(sym, isd_a[sym], scaled_adh_e[sym], a_trend[sym], imarker[n], icolour[n])
	a, d_eq, r_eq = a_trend[sym]
	trend_file.write("# Scaled_E_Adh (eV . n_interface_cluster_atoms^-1)\n#  A\tMorse interpolation:"
					 " d_eq * (exp(-2 * a * (A - r_eq)) - 2 * exp(-a * (A - r_eq)))\n")
	trend_file.write("{}\ta={:<5.5f}\td_eq={:<5.5f}\tr_eq={:<5.5f}\tR\u00B2={:<1.2f}\t\u03C4\u2264{:<1.1f}%\n"
					 .format(sym, a, d_eq, r_eq, a_r[sym], np.abs(max_deviation)))
plt.plot([z_min, 0], [z_min, 0], "k-", lw=1.5)
Display("$E_{Adh}$ (eV)", "Predicted $E_{Adh}$ (eV)", [z_min, 0], [z_min, 0], "")
trend_file.write("\n")
trend_file.close()
#-------------------------------------  Distances B ------------------------
for n, sym in enumerate(symbol):
	trend_label, b_trend[sym], b_r[sym] = trend_morse(isd_b[sym], scaled_adh_e[sym], sym, [y_min, y_max], icolour[n], imarker[n], iliner[n+1])
#y_min = min([b_trend[sym][2] for sym in symbol])*0.9
#z_min = -1*min([b_trend[sym][1] for sym in symbol]) - np.abs(min([b_trend[sym][1] for sym in symbol]))*0.1
Display("isd_b $(\\AA)$", "$E_{Adh}^{Scaled}$ $(eV \cdot atom^{-1})$", [y_min, y_max], [z_min, 0], trend_label)
#--------------------------------------- Validation ---------------------------------------
trend_file = open("Interpolation_EAdh_Trends.txt", 'a+')
for n, sym in enumerate(symbol):
	max_deviation = Validation(sym, isd_b[sym], scaled_adh_e[sym], b_trend[sym], imarker[n], icolour[n])
	b, d_eq, r_eq = b_trend[sym]
	trend_file.write("# Scaled_E_Adh (eV . n_interface_cluster_atoms^-1)\n#  B\tMorse interpolation:"
					 " d_eq * (exp(-2 * b * (B - r_eq)) - 2 * exp(-b * (B - r_eq)))\n")
	trend_file.write("{}\tb={:<5.5f}\td_eq={:<5.5f}\tr_eq={:<5.5f}\tR\u00B2={:<1.2f}\t\u03C4\u2264{:<1.1f}%\n"
					 .format(sym, b, d_eq, r_eq, b_r[sym], np.abs(max_deviation)))
plt.plot([z_min, 0], [z_min, 0], "k-", lw=1.5)
Display("$E_{Adh}$ (eV)", "Predicted $E_{Adh}$ (eV)", [z_min, 0], [z_min, 0], "")
trend_file.write("\n")
trend_file.close()