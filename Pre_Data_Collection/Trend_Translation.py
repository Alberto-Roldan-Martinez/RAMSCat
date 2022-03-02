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
	temp_energy = [float(data[i][7]) for i in range(len(data)) if float(data[i][7]) < 0]		# contains the system's energy
	i_distance = []						# contains the atom's distance to the closest neighbour
	i_atom = 0							# contains the atom index of interest
	i_coord = 0							# contains the atom's coordination
	i_gcn = 0 							# contains the atom's site generalised coordination number (the closest point)
	for i in range(len(data)):
		if float(data[i][7]) <= 0:
			if temp_energy[i] == min(temp_energy):
				i_atom = int(data[i][0])
				i_coord = int(float(data[i][1]))
				i_gcn = float(data[i][2])
				d_e_min = float(data[i][3])
				xyz_reference = np.array([float(data[i][4]), float(data[i][5]), float(data[i][6])])
	for i in range(len(data)):
		if float(data[i][7]) <= 0:
			if int(float(data[i][1])) == i_coord:
				i_distance.append(float(data[i][3]))
			else:
				vector = np.array([float(data[i][4]), float(data[i][5]), float(data[i][6])]) - xyz_reference
				d = np.sqrt(vector.dot(vector))
				i_distance.append(d_e_min + d)
	e_reference = [temp_energy[i] for i in range(len(temp_energy)) if i_distance[i] == max(i_distance)][0]
	e_coh = [i-e_reference for i in temp_energy]
#	e_coh = [i for i in t_energy]

	return i_atom, i_coord, i_gcn, i_distance, e_coh, e_reference


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
	ax.scatter3D([x0[i] for i in range(len(x0)) if xlim[0] <= x0[i] <= xlim[1] and ylim[0] <= y0[i] <= ylim[1]],
				 [y0[i] for i in range(len(y0)) if xlim[0] <= x0[i] <= xlim[1] and ylim[0] <= y0[i] <= ylim[1]],
				 [z0[i] for i in range(len(z0)) if xlim[0] <= x0[i] <= xlim[1] and ylim[0] <= y0[i] <= ylim[1]],
				 s=5, c='k', marker='o', label=trend_label)

	grid = 50
	surf_x = np.linspace(xlim[0], xlim[1], grid)
	surf_y = np.linspace(ylim[0], ylim[1], grid)
	x, y = np.meshgrid(surf_x, surf_y)
	if len(set(y0)) > 1:
		z = morse_3D([x, y], *popt)
	else:
		popt = [popt[4], popt[0], popt[1], popt[3]]
		z = morse(x, *popt)

	e_min = min([min(z[i]) for i in range(len(z))])

# masking the data beyond zmax
	z_mask_max = ma.masked_greater_equal(z, 0.0, copy=True)
	z = z_mask_max.filled(fill_value=0.0)
	z_mask_min = ma.masked_less(z, min(z0)*1.2, copy=True)
	z = z_mask_min.filled(fill_value=min(z0))

	surface = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='viridis', alpha=0.8, vmin=zlim[0], vmax=0)
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


def Extract_numeric_data(label_0, label_1, x, y, max_error):
	data_out_name = "./deviation_numeric_data"
	data_out = open(data_out_name + ".dat", "a+")
	data_out.write("# system: {} \tcoordination: {} \tmax_deviation: {:<.1} eV\n" .format(label_0, label_1, max_error))
	data_out.write("# $E_{DFT}$ $(eV \cdot atom^{\minus 1})$\t $Predicted$ $E$ $(eV \cdot atom^{\minus 1})$\n")
	for i in range(len(x)):
		data_out.write(" {:<5.5f}\t{:<5.5f}\n" .format(x[i], y[i]))
	data_out.write("\n")
	data_out.close()

def lineal(x, a, b):
	return a*x + b


def logarithm(x, a, b, c, d):
	return d/np.log(a/(a+b)) * (np.log(a) - c * np.log(a+x))


def cubic(x, a, b, c, d):
	return a*x**3 + b*x**2 + c*x + d


def morse(x, r_eq, a1, a2, d_eq):
	return d_eq * (np.exp(-2*a1*(x - r_eq)) - 2 * np.exp(-a2*(x - r_eq)))     # MORSE potential


def generalised_morse(x, r_eq, a1, a2, d_eq, m):
	return d_eq/(2*m) * ((2*m - 1) * np.exp(-2*a1*(x/r_eq - 1)) - 2*m * np.exp(-a2*(x/r_eq - 1)))     # Generalised MORSE potential: https://doi.org/10.3390/en13133323

def morse_3D(x, a1, a2, a3, a_d_eq, a_r_eq, b1, b2, b3, b_d_eq, b_r_eq):
	return a_d_eq * (np.exp(-2*a1*(x[0] - a_r_eq)) - 2 * np.exp(-(a2*(x[0] - (a_r_eq + a3/x[1]))))) +\
		   b_d_eq * (np.exp(-2*b1*(x[0] - b_r_eq))- 2 * np.exp(-(b2*(x[0] - (b_r_eq + b3/x[1])))))					# MORSE potential

def lennard_jones(x, r_eq, a, d_eq):
	return d_eq * ((r_eq/x)**(2*a) - 2*(r_eq/x)**a)     # Lennard-Jones potential


def trend_morse(x, y, symbol, xlim, colour, marker, line):
	weights = []
	r_min = min(x)
	e_min = np.abs(min(y))
	for i in range(len(x)):
#		if x[i] == max(x):
#			weights.append(0.1)
		if y[i] == min(y):
			weights.append(0.1)
			r_min = x[i]
			e_min = np.abs(y[i])
		else:
			weights.append(1)

	popt, pcov = curve_fit(morse, x, y, bounds=([r_min*0.8, 0., 0., e_min*0.8], [r_min*1.1, 20., 20., e_min*1.2]))#, sigma=weights)
	r2 = 1-np.sqrt(sum([(y[i] - morse(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
	c, a1, a2, b = popt
	trend_label = " Morse: {:.2f} * (exp(-2*{:.2f}*(d - {:.2f})) - 2 * exp(-{:.2f}*(d - {:.2f})))".format(b, a1, c, a2, c)
	print(trend_label, r2)
	x_line = np.linspace(xlim[0], xlim[1], 1500)
	y_line = morse(np.linspace(xlim[0], xlim[1], 1500), *popt)
	plt.plot(x_line, y_line, color=colour, linestyle=line, label=str(symbol) + "$\cdot R^{2}$= "+str(round(r2, 2)))
#	plt.plot(x_line, y_line, color=colour, linestyle=line, label="$R^{2}$= "+str(round(r2, 2)))
	plt.plot(x, y, marker=marker, color=colour, markersize=3, linestyle="None")
#	for i in range(len(x)):	plt.annotate(str(round(x[i], 4)), xy=(x[i]+0.001, y[i]), xytext=(x[i], y[i]))
	minima = [[x_line[i], y_line[i]] for i in range(len(y_line)) if y_line[i] == min(y_line)][0]

	return trend_label, popt, r2, minima


def trend_morse_3D(x, y, z):
	for i in range(len(x)):
		if z[i] == min(z):
			r = x[i]
	d = np.abs(min(z))

	if len(set(y)) > 1:
#			       a1, a2,   a3,  a_d_eq, a_r_eq, b1, b2,  b3,   b_d_eq, b_r_eq
		limits = ([0., 0., -r*1.2, d*0.8, r*0.8, 0., 0., -r*1.2, d*0.8, r*0.8],
				  [50, 50,  r*1.2,  d*2,  r*1.2, 50, 50,  r*1.2,  d*2,  r*1.2])
		popt, pcov = curve_fit(morse_3D, [x, y], z, bounds=limits)

		r2 = 1-np.sqrt(sum([(z[i] - morse_3D([x[i], y[i]], *popt))**2 for i in range(len(z))])/sum(i*i for i in z))
	else:
		limits = ([r, 0., 0., d], [r*2, 20., 20., d*2])
		popt, pcov = curve_fit(morse, x, z, bounds=limits)
		r2 = 1-np.sqrt(sum([(z[i] - morse(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in z))
		r_eq, a1, a2, d_eq = popt
		popt = [a1, a2, 0., d_eq, r_eq, 0., 0., 0., 0., 0.]

	return popt, r2


def Validation_3D(ele, i_coord, x0, y0, z0, popt, imarker, icolour):
	x = z0
	if popt[-2] > 0.:
		y = morse_3D(np.array([x0, y0]), *popt)
	else:
		popt = [popt[4], popt[0], popt[1], popt[3]]
		y = morse(np.array(x0), *popt)

	max_deviation = max([np.abs(y[i] - x[i]) for i in range(len(x))])
	plt.plot(x, y,  marker=imarker, color=icolour, linestyle="None", markersize=3,
			 label=str(ele) + "$^{c="+str(i_coord)+"} \cdot \\tau \leq$ " + str(round(max_deviation, 1)) + " eV")
	return max_deviation, x, list(y)

########################################################################################################################
symbol = []
i_atoms = {}
i_coords = {}
i_gcns = {}
i_distances = {}
e_coh = {}
e_reference = {}

trend_2D = {}
r2_2D = {}
trend_3D = {}
r2_3D = {}
stand_dev = {}
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol.append(str("$" + str(data[0][-2]) + "$"))					# contains the list of systems' name
	i_atoms[symbol[-1]], i_coords[symbol[-1]], i_gcns[symbol[-1]], i_distances[symbol[-1]], e_coh[symbol[-1]],\
		e_reference[symbol[-1]] = get_data(data)

x_limits = [min([min(i_distances[sym]) for sym in symbol])*0.9, 6] #max([max(i_distances[sym]) for sym in symbol])*0.8]
if max([i_gcns[sym] for sym in symbol])*1.1 <= 12.:
	y_limits = [min([i_gcns[sym] for sym in symbol])*0.9, max([i_gcns[sym] for sym in symbol])*1.1]
else:
	y_limits = [min([i_gcns[sym] for sym in symbol])*0.9, 12]
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

	label = str("c=" + str(i_coords[sym]) + " & gc=" + str(round(i_gcns[sym], 2)))
#	r_eq, a, d_eq
	trend_label, trend_2D[sym], r2_2D[sym], minima = trend_morse(i_distances[sym], e_coh[sym], label, x_limits,
														 icolour[n_colour], imarker[n_marker], iliner[n_liner])
	e_min.append(minima)
#	trend_label_2D = str(sym) + "$\cdot R^{2}$= "+"{:<1.2f}".format(float(r2_2D[sym]))
	title_label.append(str("c=" + str(i_coords[sym]) + " & gc=" + str(i_gcns[sym])))
plt.plot(x_limits, [0, 0], "k:")
# Add comments around minima
if len(e_min) > 1:
	x_position = 4
	plt.plot([e_min[0][0]-0.05, x_position*1.01], [e_min[0][1], e_min[0][1]], "k--", lw=0.5)
	plt.plot([e_min[0][0], e_min[0][0]], [e_min[0][1]-0.05, -0.1], "k--", lw=0.5)
	for i in range(1, len(e_min)):
		x0, y0 = e_min[i-1]
		x, y = e_min[i]
		plt.plot([x-0.05, x_position*1.01], [y, y], "k--", lw=0.5)
		plt.annotate("", xy=(x_position, y0), xytext=(x_position, y),
					 arrowprops=dict(arrowstyle="<->", color="k", lw=0.5))
		plt.plot([x, x], [y-0.05, -0.1], "k--", lw=0.5)
		plt.annotate("", xy=(x0, -0.5), xytext=(x, -0.5), arrowprops=dict(arrowstyle="<->", color="k", lw=0.5))
#Display("$distance$ $(\\AA)$", "$E_{Coh}^{c_{i}}$ $(eV \cdot atom^{\minus 1})$", x_limits, z_limits, "")
Display("$distance$ $(\\AA)$", "$E$ $(eV \cdot atom^{\minus 1})$", x_limits, z_limits, "")
# ------------------------------------------- 3D Display ------------------------
distances = {}
gcns = {}
coh = {}
e_mins = []
n_points = 50
for n, sym in enumerate(symbol):
	e_mins.append(min(e_coh[sym]))
	if str(i_coords[sym]) not in distances:
# Swap commented to use the 2D trends as entering points for the 3D; each with n_points
		distances[str(i_coords[sym])] = i_distances[sym]
		gcns[str(i_coords[sym])] = [i_gcns[sym] for i in range(len(i_distances[sym]))]
		coh[str(i_coords[sym])] = e_coh[sym]
#		distances[str(i_coords[sym])] = list(np.linspace(min(i_distances[sym]), max(i_distances[sym]), n_points))
#		gcns[str(i_coords[sym])] = [i_gcns[sym] for i in range(n_points)]
#		coh[str(i_coords[sym])] = list(morse(distances[str(i_coords[sym])], *trend_2D[sym]))
	else:
# Swap commented to use the 2D trends as entering points for the 3D; each with n_points
		distances[str(i_coords[sym])] += i_distances[sym]
		gcns[str(i_coords[sym])] += [i_gcns[sym] for i in range(len(i_distances[sym]))]
		coh[str(i_coords[sym])] += e_coh[sym]
#		distances[str(i_coords[sym])] += list(np.linspace(min(i_distances[sym]), max(i_distances[sym]), n_points))
#		gcns[str(i_coords[sym])] += [i_gcns[sym] for i in range(n_points)]
#		coh[str(i_coords[sym])] += list(morse(np.linspace(min(i_distances[sym]), max(i_distances[sym]), n_points), *trend_2D[sym]))
for n, coord in enumerate(distances):
	trend_3D[coord], r2_3D[coord] = trend_morse_3D(distances[str(coord)], gcns[coord], coh[coord])
	trend_label_3D = "c=" + str(coord) + "$\cdot R^{2}$= "+"{:<1.2f}".format(float(r2_3D[coord]))
	e_min = Display3D(distances[coord], gcns[coord], coh[coord], trend_3D[coord],
			  "$distance$ $(\\AA)$", "$gc$", "$E_{Coh}^{c="+coord+"}$ $(eV \cdot atom^{\minus 1})$",
					  x_limits, y_limits, z_limits, trend_label_3D)
#	e_min = Display3D(distances[coord], gcns[coord], coh[coord], trend_3D[coord],
#			  "$distance$ $(\\AA)$", "$gc$", "$E^{c="+coord+"}$ $(eV \cdot atom^{\minus 1})$",
#					  x_limits, y_limits, z_limits, trend_label_3D)

# ---------------------- Clean and get the data ----------------------------------------------
symbol = []
i_atoms = {}
i_coords = {}
i_gcns = {}
i_distances = {}
e_coh = {}
e_reference = {}
for n in range(1, len(sys.argv)):
	ifile = open(sys.argv[n]).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False and len(ifile[i].split()) > 0]
	symbol.append(str("$" + str(data[0][-2]) + "$"))					# contains the list of systems' name
	i_atoms[symbol[-1]], i_coords[symbol[-1]], i_gcns[symbol[-1]], i_distances[symbol[-1]], e_coh[symbol[-1]],\
		e_reference[symbol[-1]] = get_data(data)
# --------------------------------------- Validation ---------------------------------------
trend_file = open("Interpolation_CohesionEnergy.txt", 'w+')
max_deviation = []
x = {}
y = {}
for n, sym in enumerate(symbol):
	n_marker = n
	n_colour = n
	if n >= 2*len(icolour):
		n_colour = n - 2*len(icolour)
	elif n >= len(icolour):
		n_colour = n - len(icolour)
	if n >= len(imarker):
		n_marker = n - len(imarker)
#	print("n:", n, "\tidentifier:", sym, "\tlenght:", len(e_coh[sym]))

	deviation, x[sym], y[sym] = Validation_3D(symbol[n], int(i_coords[sym]), i_distances[sym],
						  [i_gcns[sym] for i in range(len(i_distances[sym]))], e_coh[sym], trend_3D[str(i_coords[sym])],
						  imarker[n_marker], icolour[n_colour])
	max_deviation.append(deviation)
a1, a2, a3, a_d_eq, a_r_eq, b1, b2, b3, b_d_eq, b_r_eq = trend_3D[str(i_coords[sym])]
print("Trend: ", [round(i, 5) for i in trend_3D[str(i_coords[sym])]])

if b_d_eq > 0.:
	trend_file.write("# E_Coh (eV)\t\u03C4:Maximum Absolute Error\n#" #\t\u03c3: Average Standard Deviation
				 "\t3D Morse interpolation: A + B\n"
				 "  A::\td_x_eq * ((exp(-2 * a1 * (x - r_eq)) - 2 * exp(- (a2 * (x - (r_eq + a3/y)))))\n"
				 "  B::\td_y_eq * ((exp(-2 * b1 * (y - r_eq)) - 2 * exp(- (b2 * (y - (r_eq + b3/y)))))\n")
	trend_file.write("Coordination={:d}\n\tA\td_x_eq={:<5.5f}\ta1={:<5.5f}\ta2={:<5.5f}\ta3={:<5.5f}\tr_eq={:<5.5f}\n"
				   "\tB\td_y_eq={:<5.5f}\tb1={:<5.5f}\tb2={:<5.5f}\tb3={:<5.5f}\tr_eq={:<5.5f}"
				 "\tR\u00b2={:<1.2f}  \u03C4\u2264{:<1.2f} eV\n"
				 .format(int(coord), a_d_eq, a1, a2, a3, a_r_eq, b_d_eq, b1, b2, b3, b_r_eq,
						 round(float(r2_3D[coord]), 2), np.abs(max(max_deviation))))
else:
	trend_file.write("# E_Coh (eV)\t\u03C4:Maximum Absolute Error\n#" #\t\u03c3: Average Standard Deviation
				 "\tMorse interpolation:\n"
				 "  \td_x_eq * ((exp(-2 * a1 * (x - r_eq)) - 2 * exp(- abs(a2 * (x - r_eq))))\n")
	trend_file.write("Coordination= {:d}\n\tA\td_x_eq={:<5.5f}\ta1={:<5.5f}\ta2={:<5.5f}\tr_eq={:<5.5f}"
				 "\tR\u00b2={:<1.2f}  \u03C4\u2264{:<1.2f} eV\n"
				 .format(int(coord), a_d_eq, a1, a2, a_r_eq, round(float(r2_3D[coord]), 2), np.abs(max(max_deviation))))
trend_file.close()
plt.plot(z_limits, z_limits, "k-", lw=1.5)
Display("$E_{Coh}$ $(eV \cdot atom^{\minus 1})$", "Predicted $E_{Coh}$ $(eV \cdot atom^{\minus 1})$",
		z_limits, z_limits, "")
answer = str(input("Would you like to extract numeric data from the previous plot (y/n)?\n"))
if answer == "y":
	for n, sym in enumerate(symbol):
		Extract_numeric_data(sym, int(i_coords[sym]), x[sym], y[sym], max_deviation[n])
