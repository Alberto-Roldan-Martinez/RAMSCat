'''

	USAGE: ~.py input.dat
	input: columns

'''

import sys
import numpy as np
import numpy.ma as ma
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

icolour = ["k", "b", "r", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] ## n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))]


def Display_MultiAxis(labels, x_label, x, y_labels, y, y_limits):
	fig, ax1 = plt.subplots(figsize=(12, 8), clear=True)       # prepares a figure
	ax1.set_xlabel(x_label, fontsize=18)
	ax1.tick_params(axis='x', rotation=0, labelsize=16)
#	ax1.set_xlim([0, 13])
	if len(y_labels) < 1:
		print("*** No Y axis give ***")
		exit()
	else:
		ax_position = ax1.get_position()
		ax1.set_position([ax_position.x0+(len(y_labels)-1)*0.005, ax_position.y0,
						  ax_position.x1-(len(y_labels)-1)*0.115, ax_position.y1*0.9])
		ax1.set_ylabel(y_labels[0], color=icolour[0], fontsize=16)
		ax1.tick_params(axis='y', labelcolor=icolour[0], labelsize=14)
		ax1.set_ylim(y_limits[0])
		ax1.plot(x, y[0], marker=imarker[0], linestyle=iliner[0], color=icolour[0], label=y_labels[0])
		leg_lines, leg_labels = ax1.get_legend_handles_labels()
# add labels
#    	for i in range(len(x)):
#	    	ax1.text(x[i]+0.02, y1[i], str(labels[i]), color="black", fontsize=14)
		for i in range(1, len(y_labels)):
			ax2 = ax1.twinx()           # instantiate a second axes that shares the same x-axis
			ax_position = ax2.get_position()
			ax2.set_position([ax_position.x0+(len(y_labels)-1)*0.005, ax_position.y0,
						  ax_position.x1-(len(y_labels)-1)*0.115, ax_position.y1*0.9])
			ax2.spines['right'].set_position(('axes', 1 + (i-1)*0.12)) 				# Modify the arbitrary value chage the distance between right axes
			ax2.set_ylabel(y_labels[i], color=icolour[i], fontsize=16)
			ax2.tick_params(axis='y', labelcolor=icolour[i], rotation=0, labelsize=14)
			ax2.set_ylim(y_limits[i])
			ax2.plot(x, y[i], marker=imarker[i], linestyle=iliner[i], color=icolour[i], label=y_labels[i])
			leg2_lines, leg2_labels = ax2.get_legend_handles_labels()
			leg_lines += leg2_lines
			leg_labels += leg2_labels
# add labels
#            for j in range(len(x)):
#        		ax2.text(x[i], y[i][j]+0.02, str(labels[i]), color="black", fontsize=14)
	legend = ax1.legend(leg_lines, leg_labels, loc='upper left', fontsize=14)
#	fig.tight_layout()
	plt.ion()
	plt.show()
	SaveFig()


def Display3D(labels, x, y, z, z_label):
	def plane(x, a, b, c, d, e):
		return a*np.exp(-b*x[0]) + c*np.exp(-d*x[1]) + e

#	popt, pcov = curve_fit(plane, [x, y], z)#, bounds=limits)
#	r2 = 1-np.sqrt(sum([(z[i] - plane([x, y], *popt))**2 for i in range(len(y))])/sum(i*i for i in y))
#	print(popt, "\n", r2)

	figure = plt.figure(figsize=(10, 10), clear=True)		# prepares a figure
	ax = figure.add_subplot(111, projection='3d') 			#plt.axes(projection='3d')
	ax.scatter3D(x, y, z, s=5, c='k', marker='o') #, label="$R^{2}=$ " + str(round(r2, 2)))
	z_lim = [min(z), max(z)]
# add labels
	for i in range(len(x)):
		ax.text(x[i]+0.02, y[i]+0.02, z[i]+0.02, str(labels[i]), color="black", fontsize=14)
#	spline = sp.interpolate.Rbf(x1, x2, y, function='thin_plate', smooth=5, episilon=1)

#	grid = 50
#	surf_x = np.linspace(0, max(x)*1.1, grid)
#	surf_y = np.linspace(0, max(y)*1.1, grid)
#	x, y = np.meshgrid(surf_x, surf_y)
#	z = spline(x, y)
#	z = plane([x, y], *popt)

# masking the data beyond zmax
#	z_mask_max = ma.masked_greater_equal(z, z_lim[1], copy=True)
#	z = z_mask_max.filled(fill_value=z_lim[1])
#	z_mask_min = ma.masked_less(z, z_lim[0], copy=True)
#	z = z_mask_min.filled(fill_value=z_lim[0])

#	surface = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='viridis', alpha=0.7, vmin=z_lim[0], vmax=0)
#	figure.colorbar(surface, shrink=0.25, aspect=10)

	ax.set_xlabel('coordination', fontsize=16, labelpad=10)
	ax.set_ylabel('GCN', fontsize=16)
	ax.set_zlabel(z_label, fontsize=16, labelpad=10)
#	ax.text(-2, 13, 0, str(y_label), color="black", fontsize=16)
	ax.set_xlim3d(0, max(x)*1.1)
	ax.set_ylim3d(0, max(y)*1.1)
	ax.set_zlim3d(z_lim)
	ax.tick_params(axis='both', labelsize=14)
#	ax.set_xticks([])
#	ax.set_yticks([])
#	ax.set_zticks([]) #np.linspace(0, 1, 5))
#	plt.subplots_adjust(left=0.15, right=0.9, top=0.8, bottom=0.1)
#	legend = ax1.legend(bbox_to_anchor=(0.5, 1.05), loc='upper center')
	legend = ax.legend(loc="best")
	figure.tight_layout()
	plt.grid(True)
	plt.ion()
	ax.view_init(azim=45, elev=10)
#	ax.view_init(azim=-90.01, elev=89.990)
	plt.show()
	SaveFig()
	plt.clf()

def EnergyLevels(labels, x_label, x, y_label, y, y_limit):
	inter_colum_d = 0.5
	fig, ax1 = plt.subplots(figsize=(5 + len(list(set(x))), 8), clear=True)       # prepares a figure
	ax1.set_xlabel(x_label, fontsize=16)
	print(len(list(set(x)))*inter_colum_d*3/2)
	ax1.set_xticks([inter_colum_d*3/2])	#np.arange(len(list(set(x))))
	ax1.set_xticklabels([2]) #list(set(x)))
	ax1.tick_params(axis='x', rotation=0, labelsize=14)

	x_min = min(x)
	n = 0
	for i in range(len(x)):
		if x[i] == x_min:
			if y[i] == min(y):
				ax1.plot([n + inter_colum_d, n + 1], [y[i] - min(y), y[i] - min(y)], 'b', lw=1.5)
			else:
				ax1.plot([n + inter_colum_d, n + 1], [y[i] - min(y), y[i] - min(y)], 'k', lw=0.5)

			ax1.text(n + inter_colum_d*3/2, y_limit[1] - 0.05, str(len(x)),
					 color="black", fontsize=14, ha="center", va="center",
					 bbox=dict(boxstyle="round", ec=(0.5, 0.5, 0.5), fc=(0.8, 0.8, 0.8), ))
		else:
			x_min = x[i]
			n += 1
			ax1.plot([n + 0.5, n + 1], [y[i] - min(y), y[i] - min(y)], 'k')

	ax1.set_ylabel(y_label, fontsize=16)
	ax1.tick_params(axis='y', labelsize=14)
	ax1.set_ylim(y_limit)
	fig.tight_layout()
	plt.ion()
	plt.show()
	SaveFig()



def SaveFig():
	answer = str(input("Would you like to save the figure (y/n)?\n"))
	if answer == "y":
		figure_out_name = str(input("What would it be the figure name (a word & no format)?\n"))
		plt.savefig(figure_out_name + ".svg", figsize=(12, 10), clear=True,
					bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)

######################################################################################################

ifile = open(sys.argv[1]).readlines()
data = [ifile[i].split() for i in range(len(ifile)) if len(ifile[i].split()) >= 1 and ifile[i].startswith("#") is False]
labels = [data[i][0] for i in range(len(data))]
x = [data[i][0] for i in range(len(data))]												# first column
y = [[float(data[i][j]) for i in range(len(x))] for j in range(1, len(data[0]))]				# rest of columns

y_limits = [[min(y[0])+min(y[0])*0.0005, max(y[0])-max(y[0])*0.0005],			# negative values
			[min(y[1])*0.9, max(y[1])*1.1],
			[min(y[2])*0.9, max(y[2])*1.2],
			[min(y[3])*0.99, max(y[3])*1.01]
			]
#print(labels, x, y)
#exit()

Display_MultiAxis(labels, 'Pool Size', labels, ["$E^{min}$ $(eV)$", "Time $(h)$", "Cycles $(x10^{3})$", "CPU (%)"], y, y_limits)
#Display_2axis(label, 'GCN', x2, y1, y2)

#Display3D(labels, x1, x2, y1, '$E_{eq}$ $(eV \cdot atom^{\minus 1})$')
#Display3D(labels, x, y[0], y[1], '$r_{eq}$ $(\AA)$')

#EnergyLevels(labels, 'n', x, "$E - E^{min}$ $(eV)$", y[16], [-0.05, 1])