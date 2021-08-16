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
	isd_a = []								# contains the average of the shortest distances from the  interfacial atoms in the cluster to the preferable surface site[0]
	isd_b = []								# contains the average of the shortest distances from the  interfacial atoms in the cluster to the preferable surface site[1]
	adh_e = []								# contains the DFT calculated adhesion energies
	scaled_adh_e = []
	for i in range(len(data)):
#		if float(data[i][0]) > 0:
		if float(data[i][5]) < 0:
			ic.append(float(data[i][0]))
			icc.append(float(data[i][1]))
			id.append(float(data[i][2]))
			isd_a.append(float(data[i][3]))
			isd_b.append(float(data[i][4]))
			adh_e.append(float(data[i][5]))
#				scaled_adh_e.append(float(data[i][5])/float(data[i][0])) # * float(data[i][1])))

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


def Display3D(x0, y0, z0, popt, xlabel, ylabel, zlabel, xlim, ylim, zlim, trend_label):
	figure = plt.figure(figsize=(12, 16), clear=True)		# prepares a figure
	ax = figure.add_subplot(111, projection='3d') 			#plt.axes(projection='3d')
	ax.scatter3D(x0, y0, z0, s=7, c='k', marker='o', label=trend_label)

	grid = 30
	surf_x = np.linspace(xlim[0], xlim[1], grid)
	surf_y = np.linspace(ylim[0], ylim[1], grid)
	x, y = np.meshgrid(surf_x, surf_y)
	z = morse_3D([x, y], *popt)
#	z = -np.sin(morse(x, *popt[:3]) * morse(y, *popt[-3:])) * 2   # NO!
#	z = []
#	for i in range(len(x)):
#		for j in range(len(x[i])):
#			if x[i][j] <= y[i][j]:
#				z.append(morse(x[i][j], *popt[:3]))			# >> independent Morse
#			else:
#				z.append(morse(y[i][j], *popt[-3:]))
#	z = np.reshape(z, (len(x), len(x[0])))

	surface = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='viridis', alpha=0.4)
	figure.colorbar(surface, shrink=0.25, aspect=10)
#	cset = ax.contour(x, y, z, zdir='z', offset=zlim[0], cmap='viridis')
#	cset = ax.contour(x, y, z, zdir='x', offset=max(x[-1]), cmap='viridis')
#	cset = ax.contour(x, y, z, zdir='y', offset=max(y[-1]), cmap='viridis')

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


def morse_3D(x, a, d_eq, r_eq, b, y_d_eq, y_r_eq):
#	return d_eq * (np.exp(-2*a*(x[0] - r_eq)) - 2 * np.exp(-a*(x[0] - r_eq))) +\
#		   y_d_eq * (np.exp(-2*b*(x[1] - y_r_eq)) - 2 * np.exp(-b*(x[1] - y_r_eq))) - y_d_eq		# MORSE potential
	return d_eq * (np.exp(-2*a*(x[0] - r_eq)) - 2 * np.exp(-a*(x[0] - r_eq*np.sin(x[1]/x[0])))) +\
		   y_d_eq * (np.exp(-2*b*(x[1] - y_r_eq)) - 2 * np.exp(-b*(x[1] - y_r_eq*np.sin(x[1]/x[0]))))		# MORSE potentia

def trend_morse_3D(x, y, z):
#	popt = np.array([2.03630, 1.32612, 2.14405, 1.34966, 0.31027, 3.04489])		# 2
#	popt = np.array([1.76146, 1.22805, 2.13049, 1.40225, 0.36316, 2.91604])		# 3
#	popt = np.array([1.75537, 1.34486, 2.08543, 1.10007, 0.48614, 2.95995]) 	# 4
	popt = np.array([1.5207144925510097, 1.7471659255686542, 2.224056459973962, 1.940275460489518, 0.7424901030771229, 2.4934589327284975])	#5
	r2 = 1-np.sqrt(sum([(z[i] - morse_3D([x[i], y[i]], *popt))**2 for i in range(len(x))])/sum(i*i for i in z))

	return popt, r2


def Validation_3D(ele, a, x0, y0, z0, popt, imarker, icolour):
	x = z0
	y = a * morse_3D(np.array([x0, y0]), *popt) - np.abs(max(z0)) 				        		# in eV
	max_deviation = max([(y[i] - x[i]) for i in range(len(x))])

# Add label to each point
	for i in range(len(x)):
		plt.text(x[i]+0.02, y[i]+0.02, str(i+1))

	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None",
			 label=str(ele) + "$\cdot \\tau \leq$ " + str(np.abs(round(max_deviation, 1))))

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
trend_3D = {}
r_3D = {}
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol.append(data[0][-1].split("/")[2]) # [0])					# contains the list of systems' name
	ic[symbol[-1]], icc[symbol[-1]], id[symbol[-1]], isd_a[symbol[-1]], isd_b[symbol[-1]], adh_e[symbol[-1]], scaled_adh_e[symbol[-1]] = get_data(data)
x_min = min([min(isd_a[sym]) for sym in symbol])
x_max = max([max(isd_a[sym]) for sym in symbol])*1.1
y_min = min([min(isd_b[sym]) for sym in symbol])
y_max = max([max(isd_b[sym]) for sym in symbol])*1.1
z_min = min([min(adh_e[sym]) for sym in symbol]) - np.abs(min([min(adh_e[sym]) for sym in symbol]))*0.1

#------------------------------------------- 3D Display ------------------------
for n, sym in enumerate(symbol):
	scaled_adh_e[sym] = adh_e[sym] + np.abs(max(adh_e[sym]))
	trend_3D[sym], r_3D[sym] = trend_morse_3D(isd_a[sym], isd_b[sym], scaled_adh_e[sym])
#	print("a={:<5.5f}, a_d_eq={:<5.5f}, a_r_eq={:<5.5f}\nb={:<5.5f}, b_d_eq={:<5.5f}, b_r_eq={:<5.5f}".format(*trend_3D[sym]))
	trend_label_3D = str(sym) + "$\cdot R^{2}$= "+str(round(r_3D[sym], 2))
	Display3D(isd_a[sym], isd_b[sym], scaled_adh_e[sym], trend_3D[sym], "isd_a ($\\AA$)",
			  "isd_b $(\\AA)$", "$E_{Adh}^{Scaled}$ $(eV \cdot atom^{-1})$", [x_min, x_max], [y_min, y_max], [z_min, 0], trend_label_3D)

#--------------------------------------- Validation ---------------------------------------
trend_file = open("Interpolation_EAdh.txt", 'w+')
for n, sym in enumerate(symbol):
	max_deviation = Validation_3D(sym, ic[sym], isd_a[sym], isd_b[sym], adh_e[sym], trend_3D[sym], imarker[n], icolour[n])
	a, a_d_eq, a_r_eq, b, b_d_eq, b_r_eq = trend_3D[sym]
	trend_file.write("# E_Adh (eV)\n#\t3D Morse interpolation: A + B\n")
	trend_file.write("{:<5.5f}\t{:<5.5f}{:<5.5f}\t{:<5.5f}\t{:<5.5f}\tR\u00b2={:<1.2f}\t\u03C4\u2264{:<1.1f}\n"
					 .format(a, a_d_eq, a_r_eq, b, b_d_eq, b_r_eq, round(r_3D[str(sym)], 2), np.abs(max_deviation)))
plt.plot([z_min, 0], [z_min, 0], "k-", lw=1.5)
Display("$E_{Adh}$ (eV)", "Predicted $E_{Adh}$ (eV)", [z_min, 0], [z_min, 0], "")
