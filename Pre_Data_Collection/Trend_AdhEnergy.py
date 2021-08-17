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
	popt = np.array([2.433770778924254, 1.0074054846113665, 2.18651222857168, 2.2455762462521727, 1.9020852732149833, 2.1746961810522736])		# 2
#	popt = np.array([1.9542168619948455, 1.1014439587843892, 2.3084499403864833, 1.4544180750777675, 2.5022526376198715, 1.9706719475838992])	# 3
#	popt = np.array([1.75537, 1.34486, 2.08543, 1.10007, 0.48614, 2.95995]) 	# 4
#	popt = np.array([1.9783290675048133, 1.3290659679379555, 2.2825255708618086, 1.4393564764211793, 3.2217670522828294, 1.817538002430754]		# 5

#	popt = np.array([1.9874390877442545, 1.0164703906628048, 2.2999590304359026, 1.3633612349288766, 3.6403052895520664, 1.8818775883587864])	# 7
#	popt = np.array([2.27897167607313, 1.0971623217707598, 2.231744483614781, 1.3990250677569376, 3.341951459288543, 1.8912115822136562])		# 8
#	popt = np.array([2.2342409058421655, 0.9105482462089329, 2.2610933072354733, 1.8827056640221467, 2.8573100561727145, 1.8814929505624993])	# 9
#	popt = np.array([2.36103031753729, 1.6673514647787655, 2.1110191784153174, 1.9907918148275499, 5.535263893833872, 1.522519153956438])		# 10

	r2 = 1-np.sqrt(sum([(z[i] - morse_3D([x[i], y[i]], *popt))**2 for i in range(len(x))])/sum(i*i for i in z))
	pcov = np.cov(z, morse_3D(np.array([x, y]), *popt))
	standard_deviation = sum(np.sqrt(np.diag(pcov))/np.sqrt(len(np.diag(pcov))))/len(np.diag(pcov))
	return popt, r2, standard_deviation


def Validation_3D(ele, x0, y0, z0, popt, imarker, icolour):
	x = z0
	y = morse_3D(np.array([x0, y0]), *popt) - np.abs(max(z0)) 				        		# in eV
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
stand_dev = {}
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
	trend_3D[sym], r_3D[sym], stand_dev[sym] = trend_morse_3D(isd_a[sym], isd_b[sym], scaled_adh_e[sym])
#	print("a={:<5.5f}, a_d_eq={:<5.5f}, a_r_eq={:<5.5f}\nb={:<5.5f}, b_d_eq={:<5.5f}, b_r_eq={:<5.5f}".format(*trend_3D[sym]))
	trend_label_3D = str(sym) + "$\cdot R^{2}$= "+str(round(r_3D[sym], 2))
	Display3D(isd_a[sym], isd_b[sym], scaled_adh_e[sym], trend_3D[sym], "isd_a ($\\AA$)",
			  "isd_b $(\\AA)$", "$E_{Adh}^{Scaled}$ $(eV \cdot atom^{-1})$", [x_min, x_max], [y_min, y_max], [z_min, 0], trend_label_3D)

#--------------------------------------- Validation ---------------------------------------
trend_file = open("Interpolation_EAdh.txt", 'w+')
for n, sym in enumerate(symbol):
	max_deviation = Validation_3D(sym, isd_a[sym], isd_b[sym], adh_e[sym], trend_3D[sym], imarker[n], icolour[n])
	a, a_d_eq, a_r_eq, b, b_d_eq, b_r_eq = trend_3D[sym]
	trend_file.write("# E_Adh (eV)\t\u03c3: Average Standard Deviation\t\u03C4:Maximum Absolute Error\n#"
					 "\t3D Morse interpolation: A + B\n"
					 "  A::\td_eq * (exp(-2 * a * (A - r_eq)) - 2 * exp(-a * (A - r_eq * sin(isd_b/isd_a))))\n"
					 "  B::\td_eq * (exp(-2 * b * (B - r_eq)) - 2 * exp(-b * (B - r_eq * sin(isd_b/isd_a))))\n")
	trend_file.write("{}\tA\ta={:<5.5f}\td_eq={:<5.5f} \tr_eq={:<5.5f}\n\tB\tb={:<5.5f}\td_eq={:<5.5f}\tr_eq={:<5.5f}"
					 "\tR\u00b2={:<1.2f}  \u03c3={:<1.2f}  \u03C4\u2264{:<1.2f} eV\n"
					 .format(sym, a, a_d_eq, a_r_eq, b, b_d_eq, b_r_eq, round(r_3D[sym], 2), round(stand_dev[sym], 2),
							 np.abs(max_deviation)))
plt.plot([z_min, 0], [z_min, 0], "k-", lw=1.5)
Display("$E_{Adh}$ (eV)", "Predicted $E_{Adh}$ (eV)", [z_min, 0], [z_min, 0], "")
