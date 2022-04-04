'''

	USAGE: ~.py input.dat
	input: str cc, gcn, d_eq, r_eq file with comments for labels

'''

import sys
import numpy as np
import numpy.ma as ma
import scipy as sp
import scipy.interpolate
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

icolour = ["b", "r", "k", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] ## n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))]


def Display_2axis(labels, x_label, x, y1, y2):
	fig, ax1 = plt.subplots(figsize=(8, 6), clear=True)       # prepares a figure
	ax1.set_xlabel(x_label, fontsize=16)
	ax1.tick_params(axis='x', rotation=0, labelsize=14)
	ax1.set_xlim([0, 13])
	ax1.set_ylabel('$E_{eq}$ $(eV \cdot atom^{\minus 1})$', color="blue", fontsize=16)
	ax1.tick_params(axis='y', labelcolor="blue", labelsize=14)
	ax1.set_ylim([min(y1)*1.1, max(y1)*0.9])
	ax1.plot(x, y1, "ob")
# add labels
	for i in range(len(x)):
		ax1.text(x[i]+0.02, y1[i], str(label[i]), color="black", fontsize=14)

	ax2 = ax1.twinx()           # instantiate a second axes that shares the same x-axis
	ax2.set_ylabel('$r_{eq}$ $(\AA)$', color="red", fontsize=16)
	ax2.tick_params(axis='y', labelcolor="red", rotation=0, labelsize=14)
	ax2.set_ylim([min(y2)*0.9, max(y2)*1.1])
	ax2.plot(x, y2, "sr")
# add labels
	for i in range(len(x)):
		ax2.text(x[i], y2[i]+0.02, str(label[i]), color="red", fontsize=14)

#   legend = ax1.legend(loc='best')
	fig.tight_layout()
	plt.ion()
	plt.show()
	SaveFig()


def Display3D(labels, x1, x2, y, y_label):
	def plane(x, z, a, b, c, d, e):
		return a*np.exp(-b*x[0]) + c*np.exp(-d*x[1]) + e

	popt, pcov = curve_fit(plane, [x1, x2], y)#, bounds=limits)
	r2 = 1-np.sqrt(sum([(y[i] - plane([x1[i], x2[i]], *popt))**2 for i in range(len(y))])/sum(i*i for i in y))
	print(popt, "\n", r2)

	figure = plt.figure(figsize=(10, 10), clear=True)		# prepares a figure
	ax = figure.add_subplot(111, projection='3d') 			#plt.axes(projection='3d')
	ax.scatter3D(x1, x2, y1, s=5, c='k', marker='o', label="$R^{2}=$ " + str(round(r2, 2)))
	z_lim = [min(y), max(y)]
#	spline = sp.interpolate.Rbf(x1, x2, y, function='thin_plate', smooth=5, episilon=1)

	grid = 50
	surf_x = np.linspace(0, max(x1)*1.1, grid)
	surf_y = np.linspace(0, max(x2)*1.1, grid)
	x, y = np.meshgrid(surf_x, surf_y)
#	z = spline(x, y)
	z = plane([x, y], *popt)

# masking the data beyond zmax
	z_mask_max = ma.masked_greater_equal(z, 0.0, copy=True)
	z = z_mask_max.filled(fill_value=0.0)
	z_mask_min = ma.masked_less(z, z_lim[0], copy=True)
	z = z_mask_min.filled(fill_value=z_lim[0])

	surface = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='viridis', alpha=0.7, vmin=z_lim[0], vmax=0)
#	figure.colorbar(surface, shrink=0.25, aspect=10)

	ax.set_xlabel('coordination', fontsize=16)
	ax.set_ylabel('GCN', fontsize=16)
	ax.set_zlabel(y_label, fontsize=16, labelpad=10)
	ax.set_xlim3d(0, max(x1)*1.1)
	ax.set_ylim3d(0, max(x2)*1.1)
	ax.set_zlim3d(z_lim)
	ax.tick_params(axis='both', labelsize=14)
#	ax.set_xticks([])
#	ax.set_yticks([])
#	ax.set_zticks(np.linspace(0, 1, 5))
#	plt.subplots_adjust(left=0.15, right=0.9, top=0.8, bottom=0.1)
#	legend = ax1.legend(bbox_to_anchor=(0.5, 1.05), loc='upper center')
	legend = ax.legend(loc="best")
	figure.tight_layout()
	plt.grid(True)
	plt.ion()
	ax.view_init(azim=45, elev=10)
	plt.show()
	SaveFig()
	plt.clf()


def SaveFig():
	answer = str(input("Would you like to save the figure (y/n)?\n"))
	if answer == "y":
		figure_out_name = str(input("What would it be the figure name (a word & no format)?\n"))
		plt.savefig(figure_out_name + ".svg", figsize=(12, 10), clear=True,
					bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)

######################################################################################################

ifile = open(sys.argv[1]).readlines()
data = [ifile[i].split() for i in range(len(ifile)) if ifile[i].startswith("#") is False]
label = [data[i][0] for i in range(len(data))]
x1 = [float(data[i][1]) for i in range(len(data))]
x2 = [float(data[i][2]) for i in range(len(data))]
y1 = [float(data[i][3]) for i in range(len(data))]
y2 = [float(data[i][4]) for i in range(len(data))]

#Display_2axis(label, 'coordination', x1, y1, y2)
#Display_2axis(label, 'GCN', x2, y1, y2)

Display3D(label, x1, x2, y1, '$E_{eq}$ $(eV \cdot atom^{\minus 1})$')
