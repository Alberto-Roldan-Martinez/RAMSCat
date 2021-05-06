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

def annotation(note, x0, y0, x, y):
    return plt.annotate(note, xy=(x0,y0), xycoords="data", size=14,
                               xytext=(x,y), textcoords="data",
                               arrowprops=dict(arrowstyle="<-", color="k", fc="0.6"),
                               horizontalalignment="center", verticalalignment="center")


def notes(iline):
    line = iline[-1].split("/")
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

    print(inote, len(inote))

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

    print(x, y)
    x_min = min(x)-50
    x_max = max(x)+100
    plt.plot([x_min, x_max], [x_min, x_max], "k-", lw=1.5)
    plt.plot(x, y, "rs")

    i = 0  ;   annotation(note[i], x[i], y[i], x[i]+50, y[i]-60)
    i = 1  ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]-55)
    i = 2  ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]+0)
    i = 3  ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]-75)
    i = 4  ;   annotation(note[i], x[i], y[i], x[i]+310, y[i]-35)
    i = 5  ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]+40)
    i = 6  ;   annotation(note[i], x[i], y[i], x[i]-80, y[i]+5)
    i = 7  ;   annotation(note[i], x[i], y[i], x[i]+180, y[i]-0)
    i = 8  ;   annotation(note[i], x[i], y[i], x[i]-100, y[i]+45)
#    i = 9  ;   annotation(note[i], x[i], y[i], x[i]-80, y[i]+15)
#    i = 10 ;   annotation(note[i], x[i], y[i], x[i]-70, y[i]+25)
#    i = 11 ;   annotation(note[i], x[i], y[i], x[i]-70, y[i]+65)
#    i = 12 ;   annotation(note[i], x[i], y[i], x[i]-50, y[i]+40)
#    i = 13 ;   annotation(note[i], x[i], y[i], x[i]+130, y[i]-40)
#    i = 14 ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]+17)
#    i = 15 ;   annotation(note[i], x[i], y[i], x[i]+15, y[i]+40)
#    i = 16 ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]-20)
#    i = 17 ;   annotation(note[i], x[i], y[i], x[i]+180, y[i]+20)
#    i = 18 ;   annotation(note[i], x[i], y[i], x[i]-60, y[i]+80)
#    i = 19 ;   annotation(note[i], x[i], y[i], x[i]+150, y[i]+25)
#    i = 20 ;   annotation(note[i], x[i], y[i], x[i]+250, y[i]+18)
#    i = 21 ;   annotation(note[i], x[i], y[i], x[i]+250, y[i]+40)
#    i = 22 ;   annotation(note[i], x[i], y[i], x[i]-200, y[i]-20)
#    i = 23 ;   annotation(note[i], x[i], y[i], x[i]-200, y[i]-0)
#    i = 24 ;   annotation(note[i], x[i], y[i], x[i]+10, y[i]+50)

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

    print(x, y)
    x_min = min(x)-np.abs(min(x))*0.15
    x_max = max(x)*1.15
    plt.plot([x_min, x_max], [x_min, x_max], "k-", lw=1.5)
    plt.plot(x, y, "bo")

    i = 0  ;   annotation(note[i], x[i], y[i], x[i]+0.35, y[i]+0.008)
    i = 1  ;   annotation(note[i], x[i], y[i], x[i]+0.15, y[i]-0.14)
    i = 2  ;   annotation(note[i], x[i], y[i], x[i]+0.25, y[i]-0.11)
    i = 3  ;   annotation(note[i], x[i], y[i], x[i]+0.30, y[i]+0.08)
    i = 4  ;   annotation(note[i], x[i], y[i], x[i]+0.35, y[i]+0.035)
    i = 5  ;   annotation(note[i], x[i], y[i], x[i]+0.32, y[i]-0.075)
    i = 6  ;   annotation(note[i], x[i], y[i], x[i]+0.35, y[i]-0.03)
    i = 7  ;   annotation(note[i], x[i], y[i], x[i]+0.30, y[i]-0.10)
    i = 8  ;   annotation(note[i], x[i], y[i], x[i]+0.23, y[i]-0.13)
#    i = 9  ;   annotation(note[i], x[i], y[i], x[i]+0.35, y[i]+0.07)
#    i = 10 ;   annotation(note[i], x[i], y[i], x[i]-0.14, y[i]+0.07)
#    i = 11 ;   annotation(note[i], x[i], y[i], x[i]+0.35, y[i]-0.05)
#    i = 12 ;   annotation(note[i], x[i], y[i], x[i]-0.12, y[i]+0.05)
#    i = 13 ;   annotation(note[i], x[i], y[i], x[i]-0.30, y[i])
#    i = 14 ;   annotation(note[i], x[i], y[i], x[i]-0.02, y[i]+0.12)
#    i = 15 ;   annotation(str(i) + note[i], x[i], y[i], x[i]-0.20, y[i]-0.05)
#    i = 16 ;   annotation(str(i) + note[i], x[i], y[i], x[i], y[i])
#    i = 17 ;   annotation(note[i], x[i], y[i], x[i]-0.30, y[i]-0.05)
#    i = 18 ;   annotation(note[i], x[i], y[i], x[i]-0.20, y[i]+0.00)
#    i = 19 ;   annotation(note[i], x[i], y[i], x[i]+0.18, y[i]+0.00)
#    i = 20 ;   annotation(str(i) + note[i], x[i], y[i], x[i], y[i])
#    i = 21 ;   annotation(str(i) + note[i], x[i], y[i], x[i], y[i])
#    i = 22 ;   annotation(str(i) + note[i], x[i], y[i], x[i], y[i])
#    i = 23 ;   annotation(str(i) + note[i], x[i], y[i], x[i], y[i])
#    i = 24 ;   annotation(note[i], x[i], y[i], x[i]-0.10, y[i]+0.08)

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
