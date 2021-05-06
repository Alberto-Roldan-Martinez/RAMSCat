'''

USAGE: Validate_SEnergy.py input
   input: area(m^2) SurfE(J.m^-2) and atoms coordination (i.e. from cc=3..11)

'''

import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from Library import areas, surf_energies


def Display(xlabel, ylabel, xlim, ylim, trend_label):
    plt.xlabel(str(xlabel), fontsize=14)
    plt.ylabel(str(ylabel), fontsize=14)
    plt.tick_params(axis='both', labelrotation=0, labelsize=12)               # custimise tick labels
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

def annotation(note, arrow_site, x0, y0, x, y):
    return plt.annotate(note, xy=(x0,y0), xycoords="data", size=14,
                               xytext=(x,y), textcoords="data",
                               arrowprops=dict(arrowstyle="<-", color="k", fc="0.6"),
                               horizontalalignment=arrow_site, verticalalignment="center")


def notes(iline):
    line = iline[-1].split("/")
    if line[0].startswith("#") is True:
        inote = "no"
    else:
        surface = line[0]
        size = line[1].split("_")
        if len(size) > 1:
            if "v" in size[1]:
                defect = list(size[1])[1]
                inote = str("("+surface+")$ \cdot \it p($"+size[0]+"$) \minus $"+defect)
            elif "a" in size[1]:
                defect = list(size[1])[-1]
                inote = str("("+str(surface)+")$ \cdot \it p($"+str(size[0])+"$) \plus $"+str(defect))
        else:
            inote = str("("+str(surface)+")$ \cdot \it p("+str(size[0])+")$")

    return inote

#-------------------------------------------------------------------------------------------


def Area(file_name, coord_areas):
    lines = open(file_name).readlines()

    x = []
    y = []
    note = []
    for i in range(len(lines)):
        iline = lines[i].split()
        note.append(notes(iline))
        x.append(float(iline[0]) * 1E20)
        surf_coord = [int(iline[j]) for j in range(2, len(iline)-2)]
        y.append(np.dot(coord_areas, surf_coord))                    # in Angstroms^2

    x_min = min([i for i in x + y])*0.9
    x_max = max([i for i in x + y])*1.15
    plt.plot([x_min, x_max], [x_min, x_max], "k-", lw=1.5)
    plt.plot(x, y, "rs", ms=3)
    for i in range(len(x)):
        if note[i] != "no":
            if y[i] <= x_min + (x_max - x_min)/2:
                x_margin = x_max - max([len(note[i]) for i in range(len(note))])
                y_margin = y[i]
                annotation(note[i], "left", x[i], y[i], x_margin, y_margin)
            else:
                x_margin = x_min + max([len(note[i]) for i in range(len(note))])
                y_margin = y[i]
                annotation(note[i], "right", x[i], y[i], x_margin, y_margin)
    Display("Area ($\\AA ^{2}$)", "Predicted Area ($\\AA ^{2}$)", [x_min, x_max], [x_min, x_max], "")

########################################################################

def Surf_E(file_name, coord_surf_e):
    lines = open(file_name).readlines()

    x = []
    y = []
    note = []
    for i in range(len(lines)):
        iline = lines[i].split()
        note.append(notes(iline))
        x.append(float(iline[1]))
        surf_coord = [int(iline[j]) for j in range(2, len(iline)-2)]
        y.append(np.dot(coord_surf_e, surf_coord))                    # in Angstroms^2

    x_min = min([i for i in x + y])*0.9
    x_max = max([i for i in x + y])*1.15
    plt.plot([x_min, x_max], [x_min, x_max], "k-", lw=1.5)
    plt.plot(x, y, "bo", ms=3)
    for i in range(len(x)):
        if note[i] != "no":
            if y[i] <= x_min + (x_max - x_min)/2:
                x_margin = x_max - (x_max - x_min) * max([len(note[i]) for i in range(len(note))]) / 110
                y_margin = y[i]
                annotation(note[i], "left", x[i], y[i], x_margin, y_margin)
            else:
                x_margin = x_min + (x_max - x_min) * max([len(note[i]) for i in range(len(note))]) / 110
                y_margin = y[i]
                annotation(note[i], "right", x[i], y[i], x_margin, y_margin)
    Display("$\\gamma$ ($J \cdot m^{\minus 2}$)", "Predicted $\\gamma$ ($J \cdot m^{\minus 2}$)", [x_min, x_max],
            [x_min, x_max], "")

#-----------------------------------------------------------------------------------

file_name = sys.argv[1]
element = sys.argv[2]
print("\nThe element is {}\n\n" .format(element))
coord_areas = [areas(element, i) for i in range(3, 12)]
coord_surf_e = [surf_energies(element, i) for i in range(3, 12)]

Area(file_name, coord_areas)
Surf_E(file_name, coord_surf_e)
