'''

	USAGE: ~.py input.dat
   	input: average_coordination E_Coh (eV.atom^-1) >> # Symbols Path

'''

import sys
import numpy as np
import numpy.ma as ma
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.use('TkAgg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


icolour = ["b", "r", "k", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] ## n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))]


def get_data(data):
	data.sort(key=lambda x: x[5], reverse=True)
	ic = []									# contains the number of interface_cluster_atoms
	icc = []								# contains the average coordination of the cluster atoms at the interface
	id = []									# contains the average distance from the cluster interface atoms to the surface neighbouring atoms
	isd_a = []								# contains the average of the shortest distances from the  interfacial atoms in the cluster to the preferable surface site[0]
	isd_b = []								# contains the average of the shortest distances from the  interfacial atoms in the cluster to the preferable surface site[1]
	adh_e = []								# contains the DFT calculated adhesion energies
	for i in range(len(data)):
		ic.append(float(data[i][0]))
		icc.append(float(data[i][1]))
		id.append(float(data[i][2]))
		isd_a.append(float(data[i][3]))
		isd_b.append(float(data[i][4]))
		adh_e.append(float(data[i][5]))
	if max(isd_a) >= 3:
		reference_e = [adh_e[i] for i in range(len(adh_e)) if isd_a[i] == max(isd_a)][0]
	else:
		reference_e = [adh_e[i] for i in range(len(adh_e)) if isd_b[i] == max(isd_b)][0]
	scaled_adh_e = list(np.array(adh_e) - reference_e)

	return ic, icc, id, isd_a, isd_b, adh_e, scaled_adh_e, reference_e


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
	ax.scatter3D(x0, y0, z0, s=3, c='k', marker='o', label=trend_label)

	grid = 50
	surf_x = np.linspace(xlim[0], xlim[1], grid)
	surf_y = np.linspace(ylim[0], ylim[1], grid)
	x, y = np.meshgrid(surf_x, surf_y)
	z = morse_3D([x, y], *popt)

	e_min = min([min(z[i]) for i in range(len(z))])
	for i in range(len(z)):
		for j in range(len(z[i])):
			if z[i, j] == e_min:
				x_text = x[i, j]
				y_text = y[i, j]
	ax.text(x_text, y_text, e_min*1.02, str("({:.3f},{:2.3f})".format(x_text, y_text)))

# masking the data beyond zmax
	z_mask = ma.masked_greater_equal(z, zlim[1], copy=True)
	z = z_mask.filled(fill_value=0)

	surface = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='viridis', alpha=0.4, vmin=zlim[0], vmax=0)
	figure.colorbar(surface, shrink=0.25, aspect=10)

#	cset = ax.contour(x, y, z, zdir='x', offset=max(x[-1]), cmap='viridis', alpha=0.3)
#	cset = ax.contour(x, y, z, zdir='y', offset=max(y[-1]), cmap='viridis', alpha=0.3)
	cset = ax.contour(x, y, z, zdir='z', offset=zlim[0], cmap='viridis')

	ax.set_xlabel(xlabel, fontsize=14)
	ax.set_ylabel(ylabel, fontsize=14)
	ax.set_zlabel(zlabel, fontsize=16, labelpad=10)
	ax.set_xlim3d(xlim[0], xlim[1])
	ax.set_ylim3d(ylim[0], ylim[1])
	ax.set_zlim3d(zlim[0], 0)
	ax.tick_params(axis='both', labelsize=12)
#	ax.set_xticks([])
#	ax.set_yticks([])
#	ax.set_zticks(np.linspace(0, 1, 5))
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

	return e_min

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


#def morse_3D(x, a1, a2, d_eq, r_eq, b1, b2,  y_r_eq):
#	return d_eq * ((np.exp(-2*a1*(x[0] - r_eq)) - 2 * np.exp(-a2*(x[0] - r_eq*np.sin(x[1]/x[0])))) +\
#		    (np.exp(-2*b1*(x[1] - y_r_eq)) - 2 * np.exp(-b2*(x[1] - y_r_eq*np.sin(x[1]/x[0])))))					# MORSE potential
def morse_3D(x, a1, a2, d_eq, r_eq, b1, b2, y_d_eq, y_r_eq):
	return d_eq * (np.exp(-2*a1*(x[0] - r_eq)) - 2 * np.exp(-a2*(x[0] - r_eq*np.sin(x[1]/x[0])))) +\
		   y_d_eq * (np.exp(-2*b1*(x[1] - y_r_eq)) - 2 * np.exp(-b2*(x[1] - y_r_eq*np.sin(x[1]/x[0]))))					# MORSE potentia

def trend_morse(x, y, symbol, xlim, colour, marker, line):
	popt, pcov = curve_fit(morse, x, y, bounds=([0, 0.1, 0.75], [50, 10, 5]))
	r2 = 1-np.sqrt(sum([(y[i] - morse(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
	a, b, c = popt
	trend_label = "Morse: {:.2f} * (exp(-2*{:.2f}*(id - {:.2f})) - 2 * exp(-{:.2f}*(id - {:.2f})))".format(b, a, c, a, c)
	plt.plot(np.linspace(xlim[0], xlim[1], 150), morse(np.linspace(xlim[0], xlim[1], 150), *popt), color=colour, linestyle=line)
	plt.plot(x, y, marker=marker, color=colour, markersize=3, linestyle="None",
			 label=str(symbol) + "$\cdot R^{2}$= "+str(round(r2, 2)))

	return trend_label, popt, r2


def trend_morse_3D(x, y, z):
	data_length = len(z)
	weights = [1 for i in range(len(z)-3)]
	[weights.append(0.5) for i in range(3)]
#			a1, a2, d_eq, r_eq, b1, b2, y_d_eq,  y_r_eq
	limits = ([0., 0., 0., 0.75, 0., 0., 0., 0.75], [10, 10, 10, 5, 10, 10, 50, 5])
#	initial_guess = [2., 0.8, 2., 1.5, 3., 1.8]
	popt, pcov = curve_fit(morse_3D, [x, y], z, sigma=weights, bounds=limits)

# Activate the following lines to add a 1/4 of points a z=0 around the minima forcing the function to increase faster
#	r = np.linspace(popt[4]*1.75, popt[2]*1.75, 20)
#	for i, angle in enumerate(np.arange(0, np.pi/2, np.pi/2/20)):
#		x.append(min(x)*1.5 + r[i]*np.cos(angle))
#		y.append(min(y)*1.5 + r[i]*np.sin(angle))
#		z.append(0)
#		weights.append(0.5)
#	popt, pcov = curve_fit(morse_3D, [x, y], z, p0=popt)#, sigma=weights)#, bounds=limits)

	r2 = 1-np.sqrt(sum([(z[i] - morse_3D([x[i], y[i]], *popt))**2 for i in range(data_length)])/sum(i*i for i in z))
	standard_deviation = sum(np.sqrt(np.diag(pcov)/len(np.diag(pcov))))/len(np.diag(pcov))

	return popt, r2, standard_deviation


def Validation_3D(ele, x0, y0, z0, popt, reference_e, imarker, icolour):
	x = z0
	y = morse_3D(np.array([x0, y0]), *popt) + reference_e
	max_deviation = max([np.abs(y[i] - x[i]) for i in range(len(x))])
	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None", markersize=3,
			 label=str(ele) + "$\cdot \\tau \leq$ " + str(round(max_deviation, 1)))
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
reference_adh_e = {}
trend_3D = {}
r_3D = {}
stand_dev = {}
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol.append(data[0][-1].split("/")[2]) # [0])					# contains the list of systems' name
	ic[symbol[-1]], icc[symbol[-1]], id[symbol[-1]], isd_a[symbol[-1]], isd_b[symbol[-1]], adh_e[symbol[-1]],\
	scaled_adh_e[symbol[-1]], reference_adh_e[symbol[-1]] = get_data(data)

x_min = min([min(isd_a[sym]) for sym in symbol])*0.9
x_max = max([max(isd_a[sym]) for sym in symbol])*1.1
y_min = min([min(isd_b[sym]) for sym in symbol])*0.9
y_max = max([max(isd_b[sym]) for sym in symbol])*1.1
z_min = min([min(adh_e[sym]) for sym in symbol])*1.1
z_max = max([max(adh_e[sym]) for sym in symbol])*0.9

# ------------------------------------------- 3D Display ------------------------
for n, sym in enumerate(symbol):
	trend_3D[sym], r_3D[sym], stand_dev[sym] = trend_morse_3D(isd_a[sym], isd_b[sym], scaled_adh_e[sym])
	trend_label_3D = str(sym) + "$\cdot R^{2}$= "+"{:<1.2f}".format(r_3D[sym]) #str(round(r_3D[sym], 2))

	e_min = Display3D(isd_a[sym][:len(adh_e[sym])], isd_b[sym][:len(adh_e[sym])], scaled_adh_e[sym][:len(adh_e[sym])], trend_3D[sym],
#	e_min = Display3D(isd_a[sym], isd_b[sym], scaled_adh_e[sym], trend_3D[sym],
			  "$d_{a}^{min}$ ($\\AA$)", "$d_{b}^{min}$ $(\\AA)$", "$E_{Adh}^{Scaled}$ $(eV \cdot atom^{-1})$",
			  [x_min, x_max], [y_min, y_max], [z_min, 0], trend_label_3D)
	to_predict = [round(i, 5) for i in trend_3D[sym]] + [round(reference_adh_e[sym], 5)] + [round(e_min, 5)]
	print("   To Predict:", to_predict)


# --------------------------------------- Validation ---------------------------------------
trend_file = open("Interpolation_EAdh_Sites.txt", 'w+')
for n, sym in enumerate(symbol):
	max_deviation = Validation_3D(sym, isd_a[sym][:len(adh_e[sym])], isd_b[sym][:len(adh_e[sym])], adh_e[sym],
								  trend_3D[sym], reference_adh_e[sym], imarker[n], icolour[n])
	a1, a2, a_d_eq, a_r_eq, b1, b2, b_d_eq, b_r_eq = trend_3D[sym]
#	a1, a2, d_eq, a_r_eq, b1, b2, b_r_eq = trend_3D[sym]
	trend_file.write("# E_Adh (eV)\t\u03c3: Average Standard Deviation\t\u03C4:Maximum Absolute Error\n#"
					 "\t3D Morse interpolation: A + B\n"
					 "  A::\td_eq * ((exp(-2 * a1 * (A - r_eq)) - 2 * exp(-a2 * (A - r_eq * sin(isd_b/isd_a))))\n"
					 "  B::\t\t(exp(-2 * b1 * (B - r_eq)) - 2 * exp(-b2 * (B - r_eq * sin(isd_b/isd_a)))))\n")
	trend_file.write("{}\tA\td_eq={:<5.5f}\ta1={:<5.5f}\ta2={:<5.5f}\tr_eq={:<5.5f}\n"
					   "\tB\td_eq={:<5.5f}\tb1={:<5.5f}\tb2={:<5.5f}\tr_eq={:<5.5f}"
					 "\tR\u00b2={:<1.2f}  \u03c3={:<1.2f}  \u03C4\u2264{:<1.2f} eV\n"
					 .format(sym, a_d_eq, a1, a2, a_r_eq, b_d_eq, b1, b2, b_r_eq, round(r_3D[sym], 2),
							 round(stand_dev[sym], 2), np.abs(max_deviation)))
trend_file.close()
plt.plot([z_min, 0], [z_min, 0], "k-", lw=1.5)
Display("$E_{Adh}$ (eV)", "Predicted $E_{Adh}$ (eV)", [z_min, 0], [z_min, 0], "")
