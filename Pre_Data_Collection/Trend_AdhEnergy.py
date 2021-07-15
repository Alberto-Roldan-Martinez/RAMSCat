'''

	USAGE: ~.py input.dat
   	input: average_coordination E_Coh (eV.atom^-1) >> # Symbols Path

'''

import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt



icolour = ["b", "r", "k", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] ## n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))]


def get_data(data):
	ic = []									# contains the number of interface_cluster_atoms
	icc = []								# contains the average coordination of the cluster atoms at the interface
	id = []									# contains the average distance from the cluster interface atoms to the surface neighbouring atoms
	isd = []								# contains the average of the shortest distances from the  interfacial atoms in the cluster to the preferable surface site
	adh_e = []								# contains the DFT calculated adhesion energies
	scaled_adh_e = []
	for i in range(len(data)):
		if float(data[i][0])*float(data[i][1]) > 0:
			if float(data[i][4]) < 0:
				ic.append(float(data[i][0]))
				icc.append(float(data[i][1]))
				id.append(float(data[i][2]))
				isd.append(float(data[i][3]))
				adh_e.append(float(data[i][4]))
				scaled_adh_e.append(float(data[i][4])/float(data[i][0])) # * float(data[i][1])))

	return ic, icc, id, isd, adh_e, scaled_adh_e


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


def Display3D(x, y, z, popt, xlabel, ylabel, zlabel, xlim, ylim, zlim, trend_label):
	figure = plt.figure(figsize=(12, 16), clear=True)		# prepares a figure
	ax = figure.add_subplot(111, projection='3d') 			#plt.axes(projection='3d')
	ax.scatter3D(x, y, z, s=7, c='k', marker='o')
	grid = 50
	surf_x = np.linspace(xlim[0], xlim[1], grid)
	surf_y = np.linspace(ylim[0], ylim[1], grid)
	x, y = np.meshgrid(surf_x, surf_y)
#	z = morse_3D([x, y], *popt)
	z = (morse(x, *popt[0]) + morse(y, *popt[1])) / 2
	surface = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='viridis')
	figure.colorbar(surface, shrink=0.25, aspect=10)
	cset = ax.contour(x, y, z, zdir='z', offset=zlim[0], cmap='viridis')
	cset = ax.contour(x, y, z, zdir='x', offset=max(x[-1]), cmap='viridis')
	cset = ax.contour(x, y, z, zdir='y', offset=max(y[-1]), cmap='viridis')

	ax.set_xlabel(xlabel, fontsize=14)
	ax.set_ylabel(ylabel, fontsize=14)
	ax.set_zlabel(zlabel, fontsize=14, labelpad=10)
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[0], ylim[1])
	ax.set_zlim(zlim[0], zlim[1])
	ax.tick_params(axis='both', labelsize=12)
#	plt.subplots_adjust(left=0.15, right=0.9, top=0.8, bottom=0.1)
	#legend = ax1.legend(bbox_to_anchor=(0.5, 1.05), loc='upper center')
	legend = ax.legend(loc="best")
#	figure.tight_layout()
	plt.grid(True)
	plt.ion()
	ax.view_init(azim=-135, elev=10)
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

def morse_3D(x, a, d_eq, x_r_eq, b, y_r_eq):
	return d_eq * (np.exp(-2*a*(x[0] - x_r_eq)) - 2 * np.exp(-a*(x[0] - x_r_eq))) + \
		   d_eq * (np.exp(-2*b*(x[1] - y_r_eq)) - 2 * np.exp(-b*(x[1] - y_r_eq)))			# MORSE potential

def trend_morse(x, y, symbol, xlim, colour, marker, line):
	popt, pcov = curve_fit(morse, x, y, bounds=([0, 0.1, 0.5], [50, 50, 5]))
	r2 = 1-np.sqrt(sum([(y[i] - morse(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
	a, b, c = popt
	trend_label = "Morse: {:.2f} * (exp(-2*{:.2f}*(id - {:.2f})) - 2 * exp(-{:.2f}*(id - {:.2f})))".format(b, a, c, a, c)
	plt.plot(np.linspace(xlim[0], xlim[1], 150), morse(np.linspace(xlim[0], xlim[1], 150), *popt), color=colour, linestyle=line)
	plt.plot(x, y, marker=marker, color=colour, linestyle="None", label=str(symbol) + "$\cdot R^{2}$= "+str(round(r2, 2)))

	return trend_label, popt, r2


def trend_morse_3D(x, y, z):
	popt, pcov = curve_fit(morse_3D, [x, y], z, bounds=([0, 0.1, 0.5, 0, 0.5], [50, 50, 5, 50, 5]))
	r2 = 1-np.sqrt(sum([(z[i] - morse_3D([x[i], y[i]], *popt))**2 for i in range(len(x))])/sum(i*i for i in z))

	return popt, r2


def Validation_3D(ele, x0, y0, z0, x_popt, y_popt, imarker, icolour):
	x = z0
#	y = (morse(x0, *x_popt) + morse(y0, *y_popt) / 2)          		# in eV.atom^-1
	y = morse(x0, *x_popt)
#	y = morse(y0, *y_popt)
	max_deviation = max([(y[i] - x[i])*100/x[i] for i in range(len(x))])
	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None",
			 label=str(ele) + "$\cdot \\tau \leq$ " + str(np.abs(round(max_deviation, 1))) + "%")

	return max_deviation


########################################################################################################################
symbol = []
ic = {}
icc = {}
id = {}
isd = {}
adh_e = {}
scaled_adh_e = {}
id_trend = {}
isd_trend = {}
id_r = {}
isd_r = {}
trend_3D = {}
r_3D = {}
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol.append(data[0][-1].split("/")[2]) # [0])					# contains the list of systems' name
	ic[symbol[-1]], icc[symbol[-1]], id[symbol[-1]], isd[symbol[-1]], adh_e[symbol[-1]], scaled_adh_e[symbol[-1]] = get_data(data)
x_min = min([min(id[sym]) for sym in symbol])*0.9
x_max = max([max(id[sym]) for sym in symbol])*1.3
y_min = min([min(isd[sym]) for sym in symbol])*0.9
y_max = max([max(isd[sym]) for sym in symbol])*1.3
z_min = min([min(scaled_adh_e[sym]) for sym in symbol]) - np.abs(min([min(scaled_adh_e[sym]) for sym in symbol]))*0.1
z_max = max([max(scaled_adh_e[sym]) for sym in symbol])*1.3
#------------------------------------- Interfacial Distance ------------------------
for n, sym in enumerate(symbol):
	trend_label, id_trend[sym], id_r[sym] = trend_morse(id[sym], scaled_adh_e[sym], sym, [x_min, x_max], icolour[n], imarker[n], iliner[n+1])
x_min = min([id_trend[sym][2] for sym in symbol])*0.85
z_min = -1*min([id_trend[sym][1] for sym in symbol]) - np.abs(min([id_trend[sym][1] for sym in symbol]))*0.1

Display("Interface Average Distance ($\\AA$)", "$E_{Adh}^{Scaled}$ $(eV \cdot atom^{-1})$", [x_min, x_max], [z_min, 0], trend_label)

#------------------------------------- Sorthest Distances ------------------------
for n, sym in enumerate(symbol):
	trend_label, isd_trend[sym], isd_r[sym] = trend_morse(isd[sym], scaled_adh_e[sym], sym, [y_min, y_max], icolour[n], imarker[n], iliner[n+1])
y_min = min([isd_trend[sym][2] for sym in symbol])*0.9
z_min = -1*min([isd_trend[sym][1] for sym in symbol]) - np.abs(min([isd_trend[sym][1] for sym in symbol]))*0.1
Display("isd $(\\AA)$", "$E_{Adh}^{Scaled}$ $(eV \cdot atom^{-1})$", [y_min, y_max], [z_min, 0], trend_label)

#------------------------------------------- 3D Display ------------------------
for n, sym in enumerate(symbol):
#	trend_3D[sym], r_3D[sym] = trend_morse_3D(id[sym], isd[sym], scaled_adh_e[sym])
	Display3D(id[sym], isd[sym], scaled_adh_e[sym], [id_trend[sym], isd_trend[sym]], "Interface Average Distance ($\\AA$)",
			  "isd $(\\AA)$", "$E_{Adh}^{Scaled}$ $(eV \cdot atom^{-1})$", [x_min, x_max], [y_min, y_max], [z_min, 0], sym)

#--------------------------------------- Validation ---------------------------------------
trend_file = open("Interpolation_EAdh_Trends.txt", 'w+')
for n, sym in enumerate(symbol):
	max_deviation = Validation_3D(sym, id[sym], isd[sym], adh_e[sym], id_trend[sym], isd_trend[sym], imarker[n], icolour[n])
	a, id_d_eq, id_r_eq = id_trend[sym]
	b, isd_d_eq, isd_r_eq = isd_trend[sym]
	trend_file.write("# Scaled_E_Adh (eV . n_interface_cluster_atoms^-1)\n#\tMorse interpolation:"
					 " 1/2 ( id_d_eq * (exp(-2 * a * (id - id_r_eq)) - 2 * exp(-a * (id - id_r_eq))) + "
					 "isd_d_eq * (exp(-2 * b * (isd - isd_r_eq)) - 2 * exp(-b * (isd - isd_r_eq))) )\n")
	trend_file.write("{}\ta={:<5.5f}\tid_d_eq={:<5.5f} \tid_r_eq={:<5.5f} \tR\u00B2={:<1.2f}\n"
					 "\tb={:<5.5f}\tisd_d_eq={:<5.5f}\tisd_r_eq={:<5.5f}\tR\u00B2={:<1.2f}\t\u03C4\u2264{:<1.1f}%\n"
					 .format(sym, a, id_d_eq, id_r_eq, id_r[sym], b, isd_d_eq, isd_r_eq, isd_r[sym], max_deviation))
plt.plot([z_min, 0], [z_min, 0], "k-", lw=1.5)
Display("$E_{Adh}$ (eV)", "Predicted $E_{Adh}$ (eV)", [z_min, 0], [z_min, 0], "")
