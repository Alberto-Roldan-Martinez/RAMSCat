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
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec
import matplotlib.patches as patches


icolour = ["k", "b", "r", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] ## n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))]


def Display_MultiAxis(labels, x_label, x, y_labels, y, y_limits):
	fig, ax1 = plt.subplots(figsize=(12, 8), clear=True)       # prepares a figure
	ax1.set_xlabel(x_label, fontsize=18)
	ax1.tick_params(axis='x', rotation=0, labelsize=15)
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
# 			add labels
#            for j in range(len(x)):
#        		ax2.text(x[i], y[i][j]+0.02, str(labels[i]), color="black", fontsize=14)
# add notes
#	ax2.annotate("Quasi-Newton", xy=(1, -5), xycoords="data", size=16, color='Grey', annotation_clip=False)
#	ax2.add_patch(patches.Rectangle((-0.05, -6), 3.5, 4, linewidth=1, edgecolor='cyan', facecolor='white', alpha=0.4))
#	ax2.annotate("Newtonian Dyn", xy=(4, -5), xycoords="data", size=16, color='grey', annotation_clip=False,
#				 bbox=dict(boxstyle="round, pad=0.1", fc="white", ec="white", lw=1, alpha=0.4))
#	ax2.add_patch(patches.Rectangle((3.75, -6), 1.5, 4, linewidth=1, edgecolor='cyan', facecolor='white', alpha=0.4))

	legend = ax1.legend(leg_lines, leg_labels, loc='best', fontsize=14) #upper left
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
	inter_colum_d = 0.25
	fig, ax1 = plt.subplots(figsize=(5 + len(list(set(x))), 8), clear=True)       # prepares a figure

	xtick_locator = []
	n_atoms = sorted(list(set([int(i) for i in x])))
	for n in range(len(n_atoms)):
		array = [y[i] for i in range(len(x)) if int(x[i]) == n_atoms[n]]
		for i in range(len(array)):
			if array[i] == min(array):
				ax1.plot([n_atoms[n], n_atoms[n]+1 - inter_colum_d], [array[i] - min(array), array[i] - min(array)],
						 'b', lw=1.5)
			else:
				ax1.plot([n_atoms[n], n_atoms[n]+1 - inter_colum_d], [array[i] - min(array), array[i] - min(array)],
						 'k', lw=0.5)

			ax1.text(n_atoms[n] + (1 - inter_colum_d)/2, y_limit[1] - 0.05, str(len(array)),
					 color="black", fontsize=14, ha="center", va="center",
					 bbox=dict(boxstyle="round", ec=(0.5, 0.5, 0.5), fc=(0.8, 0.8, 0.8), ))
		xtick_locator.append(n_atoms[n] + (1 - inter_colum_d)/2)

	ax1.set_xlabel(x_label, fontsize=16)
	ax1.set_ylabel(y_label, fontsize=16)
	ax1.tick_params(axis='both', rotation=0, labelsize=14)
	ax1.set_xticks(xtick_locator)
	ax1.set_xticklabels(n_atoms)
	ax1.set_ylim(y_limit)

	ax2 = ax1.twinx()           # instantiate a second axes that shares the same x-axis
	ax2.set_ylabel("$e^{ \\left( \\frac{\minus \Delta E}{k_{B} \cdot 298} \\right) }$ $(\%)$", fontsize=20, labelpad=-30)
	ax2.tick_params(axis='y', rotation=0, labelsize=14)
	ax2.set_ylim(y_limit)
	ax2.set_yticks(np.arange(0.05, 0.25, 0.05))
	y_values = [np.exp(-i/(8.617333262e-5*298))*100 for i in ax2.get_yticks()]
	ax2.set_yticklabels(["{:.1{c}}".format(i, c="e" if i < 0.1 else "f") for i in y_values])

	fig.tight_layout()
	plt.ion()
	plt.show()
	SaveFig()


def Trend(x, y, axis_min, axis_max):
	popt = np.polyfit(x, y, 1)
	poly1d_fn = np.poly1d(popt)
	r2 = 1-np.sqrt(sum([(y[i] - poly1d_fn(x[i]))**2 for i in range(len(x))])/sum(i*i for i in y))
	xx = np.linspace(axis_min, axis_max, 150)
	return popt, r2, xx, poly1d_fn(xx)


def CrossRelation(labels, x, y, x2, y2):
	fig = plt.figure(figsize=(14, 8), clear=True)
	ax = ["E_coh", "E_adh", "E_total"]
	ax[0] = fig.add_subplot(GridSpec(2, 2)[0])
	ax[1] = fig.add_subplot(GridSpec(2, 2)[2])
	ax[2] = fig.add_subplot(GridSpec(2, 2)[:, -1])
	for i in range(3):
		ax[i].tick_params(axis='both', labelsize=14)
		ax[i].axis("scaled")
		if i < 2:
			axis_max = max(x[i] + y[i])*0.9 # + x2[i] + y2[i])*0.9
			axis_min = min(x[i] + y[i])*1.1 # + x2[i] + y2[i])*1.1
		else:
			axis_max = max(x[i] + y[i])+1 # + x2[i] + y2[i])+1
			axis_min = min(x[i] + y[i])-1 # + x2[i] + y2[i])-1
		axis_step = (axis_max - axis_min)/5
		ax[i].set_xlim([axis_min, axis_max])
		ax[i].set_ylim([axis_min, axis_max])
		ax[i].xaxis.set_ticks(np.arange(axis_min, axis_max, axis_step))
		ax[i].yaxis.set_ticks(np.arange(axis_min, axis_max, axis_step))
		ax[i].xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
		ax[i].yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
		ax[i].plot([axis_min, axis_max], [axis_min, axis_max], 'k', lw=1.5)
		popt, r2, xx, yy = Trend(x[i], y[i], axis_min, axis_max)
		ax[i].plot(xx, yy)
		print("plot:", i, "equation:", popt[0], "*E_DFT +", popt[1], "R^2=", round(r2, 2))
		for j in range(len(list(set(labels)))):
			xx = [x[i][k] for k in range(len(labels)) if labels[k] == list(set(labels))[j]]
			yy = [y[i][k] for k in range(len(labels)) if labels[k] == list(set(labels))[j]]
			ax[i].plot(xx, yy, marker=imarker[j], color=icolour[j], ms=5, linestyle="None",
					   label="n= " + str(list(set(labels))[j]))
			for k in range(len(xx)):
				ax[i].annotate(k, xy=(xx[k]-0.05, yy[k]), xycoords="data")

			if len(x2) > 0:
				xx2 = [x2[i][k] for k in range(len(labels)) if labels[k] == list(set(labels))[j]]
				yy2 = [y2[i][k] for k in range(len(labels)) if labels[k] == list(set(labels))[j]]
				ax[i].plot(xx2, yy2, marker=imarker[j], color=icolour[j], ms=7, linestyle="None", markerfacecolor='none')
				for n in range(len(xx)):
					ax[i].arrow(xx[n], yy[n], xx2[n]-xx[n], yy2[n]-yy[n],# head_width=0.1, head_length=0.1,
								fc=icolour[j], ec=icolour[j])
		legend = ax[i].legend(loc="best")
	ax[0].set_xlabel("$E_{Coh}^{DFT}$ $(eV \cdot atom^{ \minus 1})$", fontsize=16)
	ax[0].set_ylabel("Predicted $E_{Coh}$ $(eV \cdot atom^{ \minus 1})$", fontsize=16)
	ax[1].set_xlabel("$E_{Adh}^{DFT}$ $(eV)$", fontsize=16)
	ax[1].set_ylabel("Predicted $E_{Adh}$ $(eV)$", fontsize=16)
	ax[2].set_xlabel("$E_{Total}^{DFT}$ $(eV)$", fontsize=16)
	ax[2].set_ylabel("Predicted $E_{Total}$ $(eV)$", fontsize=16)

	fig.tight_layout()
	plt.ion()
	plt.show()
	SaveFig()


#def Hystogram_3D(labels, x_label, x, y_label, y, y_limit):
#	figure = plt.figure(figsize=(10, 10), clear=True)		# prepares a figure
#	ax = figure.add_subplot(111, projection='3d') 			#plt.axes(projection='3d')
#	ax.scatter3D(x, y, z, s=5, c='k', marker='o') #, label="$R^{2}=$ " + str(round(r2, 2)))
#	z_lim = [min(z), max(z)]
#x, y = np.random.rand(2, 100) * 4
#hist, xedges, yedges = np.histogram2d(x, y, bins=4, range=[[0, 4], [0, 4]])
#
## Construct arrays for the anchor positions of the 16 bars.
#xpos, ypos = np.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25, indexing="ij")
#xpos = xpos.ravel()
#ypos = ypos.ravel()
#zpos = 0
#
## Construct arrays with the dimensions for the 16 bars.
#dx = dy = 0.5 * np.ones_like(zpos)
#dz = hist.ravel()
#
#ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')




def SaveFig():
	answer = str(input("Would you like to save the figure (y/n)?\n"))
	if answer == "y":
		figure_out_name = str(input("What would it be the figure name (a word & no format)?\n"))
		plt.savefig(figure_out_name + ".svg", # figsize=(12, 10), clear=True,
					bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)

def getData(fileName):
	ifile = open(fileName).readlines()
	data = [ifile[i].split() for i in range(len(ifile)) if len(ifile[i].split()) >= 1 and ifile[i].startswith("#") is False]
	labels = [data[i][0] for i in range(len(data))]
	x = [data[i][0] for i in range(len(data))]												# first column
	y = [[float(data[i][j]) for j in range(1, len(data[0]))] for i in range(len(x))]		# rest of columns

	return labels, x, y

######################################################################################################
if len(sys.argv) == 1:
	labels, x, y = getData(sys.argv[1])

#	y_limits = [[min(y[i])-np.abs(min(y[i]))*0.01, max(y[i])+np.abs(max(y[i]))*0.01] for i in range(len(y))]
	y_limits = [[min(y[0])+min(y[0])*0.0005, max(y[0])-max(y[0])*0.0005],			# negative values
				[min(y[1])*0.9, max(y[1])*1.1],
				[min(y[2])*0.9, max(y[2])*1.2],
				[min(y[3])*0.99, max(y[3])*1.01]]
	y_min = min([min(i) for i in y_limits])*1.2
	y_max = max([max(i) for i in y_limits])+0.1

#	print(labels, x, y)

	Display_MultiAxis(labels, 'GA Mutation', labels, ["$E^{min}$ $(eV)$", "Time $(h)$", "Cycles $(x10^{3})$", "CPU (%)"], y, y_limits)
#	Display_2axis(label, 'GCN', x2, y1, y2)

#	Display3D(labels, x1, x2, y1, '$E_{eq}$ $(eV \cdot atom^{\minus 1})$')
#	Display3D(labels, x, y[0], y[1], '$r_{eq}$ $(\AA)$')

#	EnergyLevels(labels, 'n', x, "$E - E^{min}$ $(eV)$", y[-1], [-0.01, 1])
else:
	labels = {}
	x = {}
	y = {}
	for n in range(1, len(sys.argv)):
		a, b, c = getData(sys.argv[n])
		labels[str(n)] = [int(a[i]) for i in range(len(a))]
		x[str(n)] = [float(i) for i in b]
		y[str(n)] = [[float(c[i][-4]) for i in range(len(c))],
					 [float(c[i][-3]) for i in range(len(c))],
					 [float(c[i][-1]) for i in range(len(c))]]					# E_Coh, E_Adh, E_Total

	if labels[str(1)] != labels[str(2)]:
		print("   The atomicity in Measured does not correspond to this in Predicted.")
		exit()

	CrossRelation(labels[str(1)], y[str(1)], y[str(2)], "", "")
#	CrossRelation(labels[str(1)], y[str(1)], y[str(3)], y[str(2)], y[str(3)])
