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
#	data.sort(key=lambda x: x[3])
	t_energy = [float(data[i][7]) for i in range(len(data)) if float(data[i][7]) < 0]		# contains the system's energy
	i_distance = []						# contains the atom's distance to the closest neighbour
	i_atom = 0							# contains the atom index of interest
	i_coord = 0							# contains the atom's coordination
	i_gcn = 0 							# contains the atom's site generalised coordination number (the closest point)
	for i in range(len(data)):
		if float(data[i][7]) < 0:
			if t_energy[i] == min(t_energy):
				i_atom = int(data[i][0])
				i_coord = int(float(data[i][1]))
				i_gcn = float(data[i][2])
				d_e_min = float(data[i][3])
				xyz_reference = np.array([float(data[i][4]), float(data[i][5]), float(data[i][6])])
	for i in range(len(data)):
		if float(data[i][7]) < 0:
			if int(float(data[i][1])) == i_coord or float(data[i][2]) == i_gcn:
				i_distance.append(float(data[i][3]))
			else:
				vector = np.array([float(data[i][4]), float(data[i][5]), float(data[i][6])]) - xyz_reference
				d = np.sqrt(vector.dot(vector))
				i_distance.append(d_e_min + d)
	e_reference = [t_energy[i] for i in range(len(t_energy)) if i_distance[i] == max(i_distance)][0]
	e_coh = [i-e_reference for i in t_energy]


	return i_atom, i_coord, i_gcn, i_distance, e_coh


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
	figure = plt.figure(figsize=(10, 10), clear=True)		# prepares a figure
	ax = figure.add_subplot(111, projection='3d') 			#plt.axes(projection='3d')
	ax.scatter3D(x0, y0, z0, s=7, c='k', marker='o', label=trend_label)

	grid = 50
	surf_x = np.linspace(xlim[0], xlim[1], grid)
	surf_y = np.linspace(ylim[0], ylim[1], grid)
	x, y = np.meshgrid(surf_x, surf_y)
	z = morse_3D([x, y], *popt)

	e_min = min([min(z[i]) for i in range(len(z))])

# masking the data beyond zmax
	z_mask = ma.masked_greater_equal(z, 0.0, copy=True)
	z = z_mask.filled(fill_value=0.0)

	surface = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='viridis', alpha=0.4, vmin=zlim[0], vmax=0)
	figure.colorbar(surface, shrink=0.25, aspect=10)
#	cset = ax.contour(x, y, z, zdir='x', offset=max(x[-1]), cmap='viridis', alpha=0.3)
#	cset = ax.contour(x, y, z, zdir='y', offset=max(y[-1]), cmap='viridis', alpha=0.3)
#	cset = ax.contour(x, y, z, zdir='z', offset=zlim[0], cmap='viridis')

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
#	legend = ax1.legend(bbox_to_anchor=(0.5, 1.05), loc='upper center')
	legend = ax.legend(loc="best")
	figure.tight_layout()
	plt.grid(True)
	plt.ion()
	ax.view_init(azim=-120, elev=10)
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


def morse(x, r_eq, a, d_eq):
	return d_eq * (np.exp(-2*a*(x - r_eq)) - 2 * np.exp(-a*(x - r_eq)))     # MORSE potential


#def morse_3D(x, a1, a2, d_eq, r_eq, b1, b2,  y_r_eq):
#	return d_eq * ((np.exp(-2*a1*(x[0] - r_eq)) - 2 * np.exp(-a2*(x[0] - r_eq*np.sin(x[1]/x[0])))) +\
#		    (np.exp(-2*b1*(x[1] - y_r_eq)) - 2 * np.exp(-b2*(x[1] - y_r_eq*np.sin(x[1]/x[0])))))					# MORSE potential
def morse_3D(x, r_eq, y_r_eq, a1, a2, d_eq, b1, b2, y_d_eq):
	return d_eq * (np.exp(-2*a1*(x[0] - r_eq)) - 2 * np.exp(-a2*(x[0] - r_eq*np.sin(x[1]/x[0])))) +\
		   y_d_eq * (np.exp(-2*b1*(x[1] - y_r_eq)) - 2 * np.exp(-b2*(x[1] - y_r_eq*np.sin(x[1]/x[0]))))					# MORSE potentia

def trend_morse(x, y, symbol, xlim, colour, marker, line):
	weights = []
	for i in range(len(x)):
		if x[i] == max(x):
			weights.append(0.1)
		else:
			weights.append(1)
	popt, pcov = curve_fit(morse, x, y, bounds=([0.75, 0, 0.], [5, 150, 100]), sigma=weights)
	r2 = 1-np.sqrt(sum([(y[i] - morse(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
	c, a, b = popt
	trend_label = "Morse: {:.2f} * (exp(-2*{:.2f}*(d - {:.2f})) - 2 * exp(-{:.2f}*(d - {:.2f})))".format(b, a, c, a, c)
	x_line = np.linspace(xlim[0], xlim[1], 150)
	y_line = morse(np.linspace(xlim[0], xlim[1], 150), *popt)
	plt.plot(x_line, y_line, color=colour, linestyle=line)
	plt.plot(x, y, marker=marker, color=colour, markersize=3, linestyle="None",
			 label=str(symbol) + "$\cdot R^{2}$= "+str(round(r2, 2)))
	minima = [[x_line[i], y_line[i]] for i in range(len(y_line)) if y_line[i] == min(y_line)][0]

	return trend_label, popt, r2, minima


def trend_morse_3D(x, y, z):
#		 r_eq, y_r_eq, a1, a2, d_eq,  b1, b2, y_d_eq
	limits = ([0.75, 0.75, 0., 0., 0., 0., 0., 0.], [5, 5, 10, 10, 50, 10, 5, 50])
	popt, pcov = curve_fit(morse_3D, [x, y], z, bounds=limits)

	r2 = 1-np.sqrt(sum([(z[i] - morse_3D([x[i], y[i]], *popt))**2 for i in range(len(z))])/sum(i*i for i in z))
	standard_deviation = sum(np.sqrt(np.diag(pcov)/len(np.diag(pcov))))/len(np.diag(pcov))

	return popt, r2, standard_deviation


def Validation_3D(ele, i_coord, x0, y0, z0, popt, imarker, icolour):
	x = z0
	y = morse_3D(np.array([x0, y0]), *popt)
	max_deviation = max([np.abs(y[i] - x[i]) for i in range(len(x))])
	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None", markersize=3,
			 label=str(ele) + "$^{c="+str(i_coord)+"} \cdot \\tau \leq$ " + str(round(max_deviation, 1)) + " eV")
	return max_deviation

########################################################################################################################
symbol = []
i_atoms = {}
i_coords = {}
i_gcns = {}
i_distances = {}
e_coh = {}

trend_2D = {}
r2_2D = {}
trend_3D = {}
r2_3D = {}
stand_dev = {}
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol.append(str("$" + str(data[0][-2]) + "$"))					# contains the list of systems' name
	i_atoms[symbol[-1]], i_coords[symbol[-1]], i_gcns[symbol[-1]], i_distances[symbol[-1]], e_coh[symbol[-1]] = get_data(data)
x_limits = [min([min(i_distances[sym]) for sym in symbol])*0.9, max([max(i_distances[sym]) for sym in symbol])*1.1]
y_limits = [min([i_gcns[sym] for sym in symbol])*0.9, max([i_gcns[sym] for sym in symbol])*1.1]
z_limits = [min([min(e_coh[sym]) for sym in symbol])*1.1, 0.1]

# ------------------------------------------- 2D Display ------------------------
title_label = []
e_min = []
for n, sym in enumerate(symbol):
	n_marker = n
	n_colour = n
	n_liner = n
	if n >= 2*len(icolour):
		n_colour = n - 2*len(icolour)
	elif n >= len(icolour):
		n_colour = n - len(icolour)
	if n >= len(imarker):
		n_marker = n - len(imarker)
	if n >= len(iliner):
		n_liner = n - len(iliner)

	label = str("cc=" + str(i_coords[sym]) + " & gcn=" + str(round(i_gcns[sym], 2)))
#	r_eq, a, d_eq
	trend_label, trend_2D[sym], r2_2D[sym], minima = trend_morse(i_distances[sym], e_coh[sym], label, x_limits,
														 icolour[n_colour], imarker[n_marker], iliner[n_liner])
	e_min.append(minima)
#	trend_label_2D = str(sym) + "$\cdot R^{2}$= "+"{:<1.2f}".format(float(r2_2D[sym]))
	title_label.append(str("cc=" + str(i_coords[sym]) + " & gcn=" + str(i_gcns[sym])))
plt.plot(x_limits, [0, 0], "k:")
# Add comments around minima
plt.plot([e_min[0][0]-0.05, 5.5], [e_min[0][1], e_min[0][1]], "k--", lw=0.5)
plt.plot([e_min[0][0], e_min[0][0]], [e_min[0][1]-0.05, -0.1], "k--", lw=0.5)
for i in range(1, len(e_min)):
	x0, y0 = e_min[i-1]
	x, y = e_min[i]
	plt.plot([x-0.05, 5.5], [y, y], "k--", lw=0.5)
	plt.annotate("", xy=(5, y0), xytext=(5, y), arrowprops=dict(arrowstyle="<->", color="k", lw=0.5))
	plt.plot([x, x], [y-0.05, -0.1], "k--", lw=0.5)
	plt.annotate("", xy=(x0, -0.5), xytext=(x, -0.5), arrowprops=dict(arrowstyle="<->", color="k", lw=0.5))
Display("$distance$ $(\\AA)$", "$E_{Coh}^{c_{i}}$ $(eV \cdot atom^{\minus 1})$", x_limits, z_limits, "")
# ------------------------------------------- 3D Display ------------------------
distances = {}
gcns = {}
coh = {}
for n, sym in enumerate(symbol):
	if str(i_coords[sym]) not in distances:
		distances[str(i_coords[sym])] = i_distances[sym]
		gcns[str(i_coords[sym])] = [i_gcns[sym] for i in range(len(i_distances[sym]))]
		coh[str(i_coords[sym])] = e_coh[sym]
	else:
		distances[str(i_coords[sym])] += i_distances[sym]
		gcns[str(i_coords[sym])] += [i_gcns[sym] for i in range(len(i_distances[sym]))]
		coh[str(i_coords[sym])] += e_coh[sym]
for n, coord in enumerate(distances):
	trend_3D[coord], r2_3D[coord], stand_dev[coord] = trend_morse_3D(distances[coord], gcns[coord], coh[coord])
	trend_label_3D = str(coord) + "$\cdot R^{2}$= "+"{:<1.2f}".format(r2_3D[coord])
	e_min = Display3D(distances[coord], gcns[coord], coh[coord], trend_3D[coord],
			  "$distance$ $(\\AA)$", "$gcn$", "$E_{Coh}^{c="+coord+"}$ $(eV \cdot atom^{\minus 1})$",
					  x_limits, y_limits, z_limits, trend_label_3D)
# --------------------------------------- Validation ---------------------------------------
trend_file = open("Interpolation_CohesionEnergy.txt", 'w+')
for n, coord in enumerate(coh):
	n_marker = n
	n_colour = n
	if n >= 2*len(icolour):
		n_colour = n - 2*len(icolour)
	elif n >= len(icolour):
		n_colour = n - len(icolour)
	if n >= len(imarker):
		n_marker = n - len(imarker)

	max_deviation = Validation_3D(coord, int(coord), distances[coord], gcns[coord], coh[coord], trend_3D[coord],
								  imarker[n], icolour[n])
	a_r_eq, b_r_eq,	a1, a2, a_d_eq, b1, b2, b_d_eq = trend_3D[coord]
#	a1, a2, d_eq, a_r_eq, b1, b2, b_r_eq = trend_3D[sym]
	trend_file.write("# E_Coh (eV)\t\u03c3: Average Standard Deviation\t\u03C4:Maximum Absolute Error\n#"
					 "\t3D Morse interpolation: A + B\n"
					 "  A::\td_eq * ((exp(-2 * a1 * (A - r_eq)) - 2 * exp(-a2 * (A - r_eq * sin(isd_b/isd_a))))\n"
					 "  B::\t\t(exp(-2 * b1 * (B - r_eq)) - 2 * exp(-b2 * (B - r_eq * sin(isd_b/isd_a)))))\n")
	trend_file.write("Coordination={:d}\n\tA\td_eq={:<5.5f}\ta1={:<5.5f}\ta2={:<5.5f}\tr_eq={:<5.5f}\n"
					   "\tB\td_eq={:<5.5f}\tb1={:<5.5f}\tb2={:<5.5f}\tr_eq={:<5.5f}"
					 "\tR\u00b2={:<1.2f}  \u03c3={:<1.2f}  \u03C4\u2264{:<1.2f} eV\n"
					 .format(int(coord), a_d_eq, a1, a2, a_r_eq, b_d_eq, b1, b2, b_r_eq, round(r2_3D[coord], 2),
							 round(stand_dev[coord], 2), np.abs(max_deviation)))
trend_file.close()
plt.plot(z_limits, z_limits, "k-", lw=1.5)
Display("$E_{Coh}$ $(eV \cdot atom^{\minus 1})$", "Predicted $E_{Coh}$ $(eV \cdot atom^{\minus 1})$", z_limits, z_limits, "")
