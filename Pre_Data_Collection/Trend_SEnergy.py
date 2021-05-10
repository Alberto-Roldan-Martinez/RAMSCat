'''

USAGE: Trend_SEnergy.py input
   input: area(m^2) SurfE(J.m^-2) and matrix of number of atoms with coordination

'''

import sys
import numpy as np
from scipy.optimize import minimize, curve_fit
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from Validate_SEnergy import Area, Surf_E


elements = []								# contains the elements forming each slab in input_file

lines = open(sys.argv[1]).readlines()
# set every array as list of each element
for i in range(len(lines)):
	iline = lines[i].split()
	elements.append(str(iline[-2]))

def get_data(lines, element):
	areas_validation = []						# contains the surface areas in m^2 that will be used for validation
	surf_e_validation = []						# contains the surface energies in J.m^2 that will be used for validation
	coord_matrix_validation = []				# contains the surface coordination
	coord_matrix_norm_validation = []			# contains the coordination matrix normalised (for surf energy)
	notes_validation = []
	areas = []									# contains the surface areas in m^2 to extract the trend
	surf_e = []									# contains the surface energies in J.m^2 to extract the trend
	coord_matrix = []
	coord_matrix_norm = []


	for i in range(len(lines)):
		iline = lines[i].split()
		if iline[-2] == element:
# data for validation plot
			if iline[-1].startswith("VAL") is True:
				areas_validation.append(float(iline[0]))				# in m^2
				surf_e_validation.append(float(iline[1]))			 	# in J . m^2
				coord_matrix_validation.append([int(iline[j]) for j in range(2, 11)])
				notes_validation.append(iline[-1][3:])
# data for interpolation and trend plot
			else:
				areas.append(float(iline[0]))				# in m^2
				surf_e.append(float(iline[1]))			 	# in J . m^2
				coord_matrix.append([int(iline[j]) for j in range(2, 11)])
				notes_validation.append("#")

# normalisation of coordinated to get the surface energy in J . m^2
	for i in range(len(coord_matrix_validation)):
		row = coord_matrix_validation[i]
		coord_matrix_norm_validation.append([float(row[j]/sum(row)) for j in range(len(row))])

	for i in range(len(coord_matrix)):
		row = coord_matrix[i]
		coord_matrix_norm.append([float(row[j]/sum(row)) for j in range(len(row))])
	coord_matrix_norm = np.reshape(coord_matrix_norm, (len(coord_matrix), len(coord_matrix)))

# solving the coordination matrix
# Areas
	coord_areas = np.linalg.solve(coord_matrix, areas)				# in m^2
#	fun = lambda x: np.linalg.norm(np.dot(coord_matrix, x) - areas)
#	sol = minimize(fun, np.zeros(len(areas)), bounds=[(0., None) for x in range(len(areas))])
#	coord_areas = sol['x']
# Surface Energies
	coord_surf_e = np.linalg.solve(coord_matrix_norm, surf_e)		# in J . m^-2
#	fun = lambda x: np.linalg.norm(np.dot(coord_matrix_norm, x) - surf_e)
#	sol = minimize(fun, np.zeros(len(surf_e)), bounds=[(0., None) for x in range(len(surf_e))])
#	coord_surf_e = sol['x']											# in J . m^-2

	return coord_areas, areas, areas_validation, coord_surf_e, surf_e, surf_e_validation,\
		   coord_matrix, coord_matrix_validation, notes_validation


def Display(xlabel, ylabel, xlim, ylim, trend_label):
	plt.xlabel(str(xlabel), fontsize=14)
	plt.ylabel(str(ylabel), fontsize=14)
	plt.tick_params(axis='both', labelrotation=0, labelsize=12)               # custimise tick labels
	plt.xticks(np.arange(int(xlim[0]), int(xlim[1]), 1))   # Xmin,Xmax,Xstep
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

def trendline_1(x, a, b, c, d):
	return a - (b + x**c)/d     		# LORENTZIAN

def trendline_2(x, a, b):
	return a*x + b                      # linear

def trend_lorentzian(x, y, element, colour, marker, line):
#          3  4  5  6  7   8    9  10 11
	weights = [1, 1, 1, 1, 1, 0.5, 0.5, 1, 1]
	popt1, pcov1 = curve_fit(trendline_1, x, y, sigma=weights)
	r2 = 1-np.sqrt(sum([(y[i] - trendline_1(x[i], *popt1))**2 for i in range(len(x))])/sum(i*i for i in y))
	a, b, c, d = popt1
	trend_label = "Lorentzian: a, b, c, d =" + str(round(a, 5)) + ", " +\
				  str(round(b, 5)) + ", " + str(round(c, 5)) + ", " + str(round(d, 5))
	plt.plot(x, y, marker=marker, color=colour, linestyle="None")
	plt.plot(np.linspace(0, 12, 150), trendline_1(np.linspace(0, 12, 150), *popt1), color=colour, linestyle=line,
			 label=str(element) + "$\cdot R^{2}$= "+str(round(r2, 2)))

	return trend_label, popt1

def trend_lineal(x, y, element, colour, marker, line):
	weights = [1, 1, 1, 1, 1, 0.5, 0.5, 1, 1]
	popt2, pcov2 = curve_fit(trendline_2, x, y, sigma=weights)#, absolute_sigma=True)#, bounds=([0, 0, 0, 20000], [50, 150, 150, 40000]))
	r2 = 1-np.sqrt(sum([(y[i] - trendline_2(x[i], *popt2))**2 for i in range(len(x))])/sum(i*i for i in y))
	a, b = popt2
	trend_label = "a*x + b =" + str(round(a, 5)) + "x +" + str(round(b, 5))
	plt.plot(x, y, marker=marker, color=colour, linestyle="None")
	plt.plot(np.linspace(0, 12, 150), trendline_2(np.linspace(0, 12, 150), *popt2), color=colour, linestyle=line,
			 label=str(element) + "$\cdot R^{2}$= "+str(round(r2, 2)))

	return trend_label, popt2

#####################################################################################

icolour = ["b", "r", "k", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] ## n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))]

# x axis array for plotting the coordination from cc=2..11
x_axis = np.arange(3, 12, 1)
for i, ele in enumerate(set(elements)):
	returned_data = get_data(lines, ele)
	coord_areas, areas, areas_validation = returned_data[0:3]
	coord_surf_e, surf_e, surf_e_validation = returned_data[3:6]
	coord_matrix, coord_matrix_validation, notes_validation = returned_data[6:]


# Area trend
	y = coord_areas * 1E20										# in Angstroms^2
	trend_label, popt_lorentzian = trend_lorentzian(x_axis, y, ele, icolour[i], imarker[i], iliner[i+1])
	trend_label, popt_lineal = trend_lineal(x_axis, y, ele, icolour[i], imarker[-i], iliner[i+1])
	areas_excel = [1.016E-19, 9.990E-20, 9.618E-20, 9.340E-20, 9.252E-20, 8.702E-20, 7.54E-20, 6.238E-20, 5.847E-20]
#	print(y, "\n", [float(i)*1E20 for i in areas_excel])
	plt.plot(x_axis, [float(i)*1E20 for i in areas_excel], "ko", ms=2)
#	Display("Coordination", "Area ($\\AA ^{2} \cdot atom^{\minus 1}$)", [0, 12.15], [min(y)-np.abs(min(y)*0.15), max(y)*1.15], trend_label)

	print(notes_validation)

	Area(ele, [j*1E20 for j in areas+areas_validation], coord_matrix+coord_matrix_validation, popt_lineal, notes_validation)

# Surface energy trend
#	y = coord_surf_e										# in J . m^-2
##	trend_label, popt_lorentzian = trend_lorentzian(x_axis, y, ele, icolour[i], imarker[i], iliner[i+1])
#	trend_label, popt_lineal = trend_lineal(x_axis, y, ele, icolour[i], imarker[i], iliner[i+1])
#	surf_e_excel = [2.796, 2.525, 2.383, 1.917, 1.785, 1.513, 1.368, 0.813, 0.374]
##	print(y, "\n", surf_e_excel)
#	plt.plot(x_axis, surf_e_excel, "ko", ms=2)
#	Display("Coordination", "$\\gamma$ ($J \cdot m^{\minus 2}$)", [0, 12.15], [min(coord_surf_e)-np.abs(min(coord_surf_e)*0.15), max(coord_surf_e)*1.15], trend_label)
#
#	toeV = 1.60218E+19
#	y = coord_surf_e * toeV * coord_areas        				# energies in eV . atom^-1
#	trend_label, popt_lineal = trend_lineal(x_axis, y, ele, icolour[i], imarker[i], iliner[i+1])
#	Display("Coordination", "$\\gamma$ ($eV \cdot atom^{\minus 1}$)", [0, 12.15], [min(y)-np.abs(min(y)*0.15), max(y)*1.15], trend_label)

