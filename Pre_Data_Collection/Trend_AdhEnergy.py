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
	idr = []								# contains the shortest distance from the reference (FIRST interfacial atom in the cluster) to the a surface neighbour
	adh_e = []								# contains the DFT calculated adhesion energies
	scaled_adh_e = []
	for i in range(len(data)):
		if float(data[i][0])*float(data[i][1]) > 0:
			if float(data[i][4]) < 0:
				ic.append(float(data[i][0]))
				icc.append(float(data[i][1]))
				id.append(float(data[i][2]))
				idr.append(float(data[i][3]))
				adh_e.append(float(data[i][4]))
				scaled_adh_e.append(float(data[i][4])/(float(data[i][0]) * float(data[i][1])))

	return ic, icc, id, idr, adh_e, scaled_adh_e


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


def Display3D(x, y, z, xlabel, ylabel, zlabel, xlim, ylim, zlim, trend_label):
	figure = plt.figure()       						# prepares a figure
	ax = figure.add_subplot(111, projection='3d') 		#plt.axes(projection='3d')
	ax.scatter3D(x, y, z, cmap='viridis', label=trend_label)

#    x = np.linspace(popt_Mg[2]*0.85, popt_Mg[2]*1.5, grid)
#    y = np.linspace(popt_O[2]*0.85, popt_O[2]*1.5, grid)
#    x,y = np.meshgrid(x,y)
#    z = function(x,popt_Mg[0],popt_Mg[1],popt_Mg[2])+function(y,popt_O[0],popt_O[1],popt_O[2])
#	 surface = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='viridis')
#    figure.colorbar(surface, shrink=0.25, aspect=10)
#    cset = ax.contour(x, y, z, zdir='z', offset=zmin*1.1, cmap='viridis')
#    cset = ax.contour(x, y, z, zdir='x', offset=max(x[-1]), cmap='viridis')
#    cset = ax.contour(x, y, z, zdir='y', offset=max(y[-1]), cmap='viridis')

	ax.set_xlabel(xlabel, fontsize=14)
	ax.set_ylabel(ylabel, fontsize=14)
	ax.set_zlabel(zlabel, fontsize=14)
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[0], ylim[1])
	ax.set_zlim(zlim[0], zlim[1])
	ax.tick_params(axis='both', labelrotation=0, labelsize=12)
#	plt.subplots_adjust(left=0.15, right=0.9, top=0.8, bottom=0.1)
	#legend = ax1.legend(bbox_to_anchor=(0.5, 1.05), loc='upper center')
	legend = ax.legend(loc="best")
	figure.tight_layout()
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

def trend_morse(x, y, symbol, colour, marker, line):
	popt, pcov = curve_fit(morse, x, y, bounds=([0, 0.1, 0.5], [50, 50, 5]))
	r2 = 1-np.sqrt(sum([(y[i] - morse(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
	a, b, c = popt
	trend_label = "Morse: {:.2f} * (exp(-2*{:.2f}*(id - {:.2f})) - 2 * exp(-{:.2f}*(id - {:.2f})))".format(b, a, c, a, c)
	plt.plot(np.linspace(0, 12, 150), morse(np.linspace(0, 12, 150), *popt), color=colour, linestyle=line)
	plt.plot(x, y, marker=marker, color=colour, linestyle="None", label=str(symbol) + "$\cdot R^{2}$= "+str(round(r2, 2)))

	return trend_label, popt, r2


def Validation(ele, x0, y0, popt, imarker, icolour):
	x = y0
	y = [morse(i, *popt) for i in x0]          		# in eV.atom^-1
	max_deviation = max([(y[i] - x[i])*100/x[i] for i in range(len(x))])
	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None",
			 label=str(ele) + "$\cdot \\tau \leq$ " + str(np.abs(round(max_deviation, 1))) + "%")

	return max_deviation


########################################################################################################################
trend = {}
r = {}
z_min = 1
z_max = 1
y_min = 1
y_max = 1
x_min = 1
x_max = 1
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol = data[0][-1].split("/")[2] # [0]					# contains the list of systems' name
	ic, icc, id, idr, adh_e, scaled_adh_e = get_data(data)
	if x_min > min(id):
		x_min = min(id)
	if x_max < max(id):
		x_max = max(id)
	if y_min > min(idr):
		y_min = min(idr)
	if y_max < max(idr):
		y_max = max(idr)
	if z_min > min(scaled_adh_e):
		z_min = min(scaled_adh_e)
	if z_max < max(scaled_adh_e):
		z_max = max(scaled_adh_e)
#------------------------------------- Interfacial Distance ------------------------
	trend_label, trend[symbol], r[symbol] = trend_morse(id, scaled_adh_e, symbol, icolour[n-1], imarker[n-1], iliner[n])
x_min = x_min*0.8
x_max = x_max*1.5
y_min = y_min*0.8
y_max = y_max*1.5
z_min = z_min - np.abs(z_min*0.1)
z_max = z_max*1.05
Display("Interface Average Distance ($\\AA$)", "$E_{Adh}^{Scaled}$ $(eV \cdot atom^{-1})$", [x_min, x_max], [z_min, z_max], trend_label)

#------------------------------------------- 3D Display ------------------------
x = []
y = []
z = []
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol = data[0][-1].split("/")[2] # [0]					# contains the list of systems' name
	ic, icc, id, idr, adh_e, scaled_adh_e = get_data(data)
	x.append(id)
	y.append(idr)
	z.append(scaled_adh_e)

Display3D(x, y, z, "Interface Average Distance ($\\AA$)", "The Shortest Interfacial Distance ($\\AA$)",
		  "$E_{Adh}^{Scaled}$ $(eV \cdot atom^{-1})$", [x_min, x_max], [y_min, y_max], [z_min, z_max], symbol)

z_min = 1
z_max = -1
trend_file = open("Interpolation_EAdh_Trends.txt", 'w+')
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol = data[0][-1].split("/")[2] #[0]					# contains the list of systems' name
	ic, icc, id, idr, adh_e, scaled_adh_e = get_data(data)
	if z_min > min(adh_e):
		z_min = min(adh_e)
	if z_max < max(adh_e):
		z_max = max(adh_e)
#--------------------------------------- Validation ---------------------------------------
	max_deviation = Validation(symbol[:2], id, adh_e, trend[symbol], imarker[n-1], icolour[n-1])
	trend_file.write("# Scaled_E_Adh (eV . n_interface_cluster_atoms^-1)\tMorse interpolation: d_eq * (np.exp(-2*a*(x - r_eq)) - 2 * np.exp(-a*(x - r_eq)))\n")
	trend_file.write("{}\ta={:<5.5f}\td_eq={:<5.5f}\tr_eq={:<5.5f}\tR\u00B2={:<1.2f}\t\u03C4\u2264{:<1.1f}%\n"
					 .format(symbol, trend[symbol][0], trend[symbol][1], trend[symbol][2], r[symbol], np.abs(round(max_deviation, 1))))
z_min = z_min - np.abs(z_min*0.1)
z_max = z_max + np.abs(z_max*0.1)
plt.plot([z_min, z_max], [z_min, z_max], "k-", lw=1.5)
Display("$E_{Adh}$ (eV)", "Predicted $E_{Adh}$ (eV)", [z_min, z_max], [z_min, z_max], "")
