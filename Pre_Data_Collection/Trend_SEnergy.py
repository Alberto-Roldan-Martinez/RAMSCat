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

data = np.loadtxt(sys.argv[1], comments="#")                    # import and reads data
x_axis = np.arange(3, 12, 1)

areas = [float(i) for i in data[:, 0]]							# in m^2
surf_e = [float(i) for i in data[:, 1]]							# in J . m^2
coord_matrix = data[:, 2:].astype(int)
coord_matrix_norm = []

for i in range(len(coord_matrix)):
	row = coord_matrix[i, :]
	coord_matrix_norm.append([float(row[j]/sum(row)) for j in range(len(row))])
coord_matrix_norm = np.reshape(coord_matrix_norm, (len(coord_matrix), len(coord_matrix)))

# solving the coordination matrix
# Areas
coord_areas = np.linalg.solve(coord_matrix, areas)				# in m^2 . atom^-1
# Surface Energies
fun = lambda x: np.linalg.norm(np.dot(coord_matrix_norm, x) - surf_e)
sol = minimize(fun, np.zeros(len(surf_e)), bounds=[(0., None) for x in range(len(surf_e))])
coord_surf_e = sol['x']											# in J . m^-2 . atom^-1


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

def trend_lorentzian(x, y):
#          3  4  5  6  7   8    9  10 11
	weights = [1, 1, 1, 1, 1, 0.5, 0.5, 1, 1]
	popt1, pcov1 = curve_fit(trendline_1, x, y, sigma=weights)
	r2 = 1-np.sqrt(sum([(y[i] - trendline_1(x[i], *popt1))**2 for i in range(len(x))])/sum(i*i for i in y))
	a, b, c, d = popt1
	trend_label = "Lorentzian: a, b, c, d =" + str(round(a, 5)) + ", " +\
				  str(round(b, 5)) + ", " + str(round(c, 5)) + ", " + str(round(d, 5))
	plt.plot(x, y, "rs")
	plt.plot(np.linspace(0, 12, 150), trendline_1(np.linspace(0, 12, 150), *popt1),
			 "r--", label="$R^{2}$= "+str(round(r2, 2)))

	return trend_label

def trend_lineal(x, y):
	weights = [1, 1, 1, 1, 1, 0.5, 0.5, 1, 1]
	popt2, pcov2 = curve_fit(trendline_2, x, y, sigma=weights)#, absolute_sigma=True)#, bounds=([0, 0, 0, 20000], [50, 150, 150, 40000]))
	r2 = 1-np.sqrt(sum([(y[i] - trendline_2(x[i], *popt2))**2 for i in range(len(x))])/sum(i*i for i in y))
	a, b = popt2
	trend_label = "a*x + b =" + str(round(a, 5)) + "x +" + str(round(b, 5))
	plt.plot(x, y, "bo")
	plt.plot(np.linspace(0, 12, 150), trendline_2(np.linspace(0, 12, 150), *popt2),
			 "b--", label="$R^{2}$= "+str(round(r2, 2)))

	return trend_label
'''

Area trend

'''
y = coord_areas * 1E20										# in Angstroms^2
trend_label = trend_lorentzian(x_axis, y)
trend_label = trend_lineal(x_axis, y)
Display("Coordination", "Area ($\\AA ^{2} \cdot atom^{\minus 1}$)", [0, 12.15], [min(y)-np.abs(min(y)*0.15), max(y)*1.15], trend_label)
'''

Surface energy trend

'''
y = coord_surf_e										# in J . m^-2
trend_label = trend_lineal(x_axis, y)
#trend_label = trend_lorentzian(x_axis, y)
Display("Coordination", "$\\gamma$ ($J \cdot m^{\minus 2}$)", [0, 12.15], [min(coord_surf_e)-np.abs(min(coord_surf_e)*0.15), max(coord_surf_e)*1.15], trend_label)

toeV = 1.60218E+19
y = coord_surf_e * toeV * coord_areas        				# energies in eV . atom^-1
trend_label = trend_lineal(x_axis, y)
Display("Coordination", "$\\gamma$ ($eV \cdot atom^{\minus 1}$)", [0, 12.15], [min(y)-np.abs(min(y)*0.15), max(y)*1.15], trend_label)

