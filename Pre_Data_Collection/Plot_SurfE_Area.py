'''

USAGE: Plot_area.py coord_area_data all_areas_data

'''



import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

data_coord = np.loadtxt(sys.argv[1], comments="#")                    # import and reads data
data = open(sys.argv[2]).readlines()
data = [data[i].split() for i in range(len(data))]

x0 = [i for i in data_coord[:, 0] if i != 12]
y0 = [data_coord[i, 1]*1E20 for i in range(len(data_coord[:, 1])) if data_coord[i, 0] != 12]         # 1= Areas in m^2
x12 = [i for i in data_coord[:, 0] if i == 12]
y12 = [data_coord[i, 1]*1E20 for i in range(len(data_coord[:, 1])) if data_coord[i, 0] == 12]

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

def trendline_1(x, a, b, c, d):
    return a - (b + x**c)/d     # LORENTZIAN

def trendline_2(x, a, b):
    return a*x + b                      # linear


# Find the Trendline
#          3  4  5  6  7   8    9  10 11
#weights = [1, 1, 1, 1, 1, 0.2, 0.1, 1, 1]
popt1, pcov1 = curve_fit(trendline_1, x0, y0) #, sigma=weights)#, absolute_sigma=True)#, bounds=([0, 0, 0, 20000], [50, 150, 150, 40000]))
r2 = 1-np.sqrt(sum([(y0[i] - trendline_1(x0[i], *popt1))**2 for i in range(len(x0))])/sum(i*i for i in y0))
a, b, c, d = popt1
trend_label = "Lorentzian: a, b, c, d =" + str(round(a, 5)) + ", " +\
                  str(round(b, 5)) + ", " + str(round(c, 5)) + ", " + str(round(d, 5))
# What to plot?
plt.plot(x0, y0, "rs")
if int(x12[0]) == 12:
    plt.plot(x12, y12, "rs", fillstyle="none")
plt.plot(np.linspace(0, 12, 150), trendline_1(np.linspace(0, 12, 150), *popt1), "r--", label="$R^{2}$= "+str(round(r2, 2)))

Display("Coordination", "Area ($\\AA ^{2} \cdot atom^{\minus 1}$)", [0, 12.15], [min(y0)-1, max(y0)*1.15], trend_label)

#-------------------------------------------------------------------------------------------

# What to plot?
coord = [float(data[i][4]) for i in range(len(data))]               # coordination to get predicted area
natoms = [float(data[i][3]) for i in range(len(data))]              # number of atoms with such coordination
x = [float(data[i][1])*1E20 for i in range(len(data))]              # real Area
y = [trendline_1(coord[i], *popt1)*natoms[i] for i in range(len(coord))]

plt.plot([-50, max(x)+100], [-50, max(x)+100], "k-", lw=1.5)
plt.plot(x, y, "rs")
note = []
for i in range(len(x)):
    line = data[i][0].split("/")
    surface = line[0]
    size = line[1].split("_")
    if len(size) > 1:
        if "v" in size[1]:
            defect = list(size[1])[1]
            note.append(str("("+surface+")$ \cdot \it p($"+size[0]+"$) \minus $"+defect))
        elif "a" in size[1]:
            defect = list(size[1])[-1]
            note.append(str("("+str(surface)+")$ \cdot \it p($"+str(size[0])+"$) \plus $"+str(defect)))
    else:
        note.append(str("("+str(surface)+")$ \cdot \it p("+str(size[0])+")$"))

i = 0  ;   annotation(note[i], x[i], y[i], x[i]+50, y[i]-60)
i = 1  ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]-55)
i = 2  ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]+0)
i = 3  ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]-75)
i = 4  ;   annotation(note[i], x[i], y[i], x[i]+310, y[i]-35)
i = 5  ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]+40)
i = 6  ;   annotation(note[i], x[i], y[i], x[i]-80, y[i]+5)
i = 7  ;   annotation(note[i], x[i], y[i], x[i]+180, y[i]-0)
i = 8  ;   annotation(note[i], x[i], y[i], x[i]-100, y[i]+45)
i = 9  ;   annotation(note[i], x[i], y[i], x[i]-80, y[i]+15)
i = 10 ;   annotation(note[i], x[i], y[i], x[i]-70, y[i]+25)
i = 11 ;   annotation(note[i], x[i], y[i], x[i]-70, y[i]+65)
i = 12 ;   annotation(note[i], x[i], y[i], x[i]-50, y[i]+40)
i = 13 ;   annotation(note[i], x[i], y[i], x[i]+130, y[i]-40)
i = 14 ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]+17)
i = 15 ;   annotation(note[i], x[i], y[i], x[i]+15, y[i]+40)
i = 16 ;   annotation(note[i], x[i], y[i], x[i]+300, y[i]-20)
i = 17 ;   annotation(note[i], x[i], y[i], x[i]+180, y[i]+20)
i = 18 ;   annotation(note[i], x[i], y[i], x[i]-60, y[i]+80)
i = 19 ;   annotation(note[i], x[i], y[i], x[i]+150, y[i]+25)
i = 20 ;   annotation(note[i], x[i], y[i], x[i]+250, y[i]+18)
i = 21 ;   annotation(note[i], x[i], y[i], x[i]+250, y[i]+40)
i = 22 ;   annotation(note[i], x[i], y[i], x[i]-200, y[i]-20)
i = 23 ;   annotation(note[i], x[i], y[i], x[i]-200, y[i]-0)
i = 24 ;   annotation(note[i], x[i], y[i], x[i]+10, y[i]+50)

Display("Area ($\\AA ^{2}$)", "Predicted Area ($\\AA ^{2}$)", [-30, max(x)*1.1], [-30, max(x)*1.1], "")

########################################################################

# What to plot?
x0 = [i for i in data_coord[:, 0] if i != 12]
y0 = [data_coord[i, 2] for i in range(len(data_coord[:, 2])) if data_coord[i, 0] != 12]         # 2= Surface energies J.m^2
x12 = [i for i in data_coord[:, 0] if i == 12]
y12 = [data_coord[i, 2] for i in range(len(data_coord[:, 2])) if data_coord[i, 0] == 12]

# Find the Trendline
popt2, pcov2 = curve_fit(trendline_1, x0, y0)# , p0=[0, 0, 0, 1]) #, bounds=([0, 0, 0, 1], [50, 150, 150, 40]))
r2 = 1-np.sqrt(sum([(y0[i] - trendline_1(x0[i], *popt2))**2 for i in range(len(x0))])/sum(i*i for i in y0))
a, b, c, d = popt2
trend_label = "Lorentzian: a, b, c, d =" + str(round(a, 5)) + ", " + str(round(b, 5)) + ", " + str(round(c, 5)) +\
              ", " + str(round(d, 5))
plt.plot(x0, y0, "bo")
if int(x12[0]) == 12:
    plt.plot(x12, y12, "bo", fillstyle="none")
plt.plot(np.linspace(0, 12, 150), trendline_1(np.linspace(0, 12, 150), *popt2), "b--", label="$R^{2}$= " + str(round(r2, 2)))

Display("Coordination", "$\\gamma$ ($J \cdot m^{\minus 2}$)", [0, 12.15], [min(y0)-min(y0)*1.1, max(y0)*1.15], trend_label)

#-------------------------------------------------------------------------------------------------
# What to plot?
toeV = 1.60218E19

x0 = [i for i in data_coord[:, 0] if i != 12]
y0 = [data_coord[i, 2]*data_coord[i, 1]*toeV for i in range(len(data_coord[:, 2])) if data_coord[i, 0] != 12]         # 2= Surface energies J.m^2
x12 = [i for i in data_coord[:, 0] if i == 12]
y12 = [data_coord[i, 2]*data_coord[i, 1]*toeV for i in range(len(data_coord[:, 2])) if data_coord[i, 0] == 12]

# Find the Trendline
#          3  4  5  6  7   8    9  10 11
#weights = [1, 1, 1, 1, 1, 0.2, 0.1, 1, 1]
popt2, pcov2 = curve_fit(trendline_2, x0, y0)   #, sigma=weights, absolute_sigma=True)#, bounds=([0, 0, 0, 20000], [50, 150, 150, 40000]))
r2 = 1-np.sqrt(sum([(y0[i] - trendline_2(x0[i], *popt2))**2 for i in range(len(x0))])/sum(i*i for i in y0))
a, b = popt2
trend_label = "a*x + b =" + str(a) + "x +" + str(b)
plt.plot(x0, y0, "bo")
if int(x12[0]) == 12:
    plt.plot(x12, y12, "bo", fillstyle="none")
plt.plot(np.linspace(0, 12, 150), trendline_2(np.linspace(0, 12, 150), *popt2), "b--", label="$R^{2}$= " + str(round(r2, 2)))

Display("Coordination", "$\\gamma$ ($eV \cdot atom^{\minus 1}$)", [0, 12.15], [-0.2, max(y0)*1.15], trend_label)

x = [float(data[i][2])*float(data[i][1])/float(data[i][3])*toeV for i in range(len(data))]
y2 = [trendline_2(coord[i], *popt2) for i in range(len(data))]
y = []
for i in range(len(data)):
    atoms = [float(data[i][j+5]) for j in range(len(y0)) if float(data[i][j+5]) != 0.]
    y.append(sum([y0[j]*float(data[i][j+5]) for j in range(len(y0)) if float(data[i][j+5]) != 0.])/natoms[i])

xx = []
yy = []

i = 0  ;   annotation(note[i], x[i], y[i], x[i]+0.35, y[i]+0.008); xx.append(x[i]); yy.append(y[i])
i = 1  ;   annotation(note[i], x[i], y[i], x[i]+0.15, y[i]-0.14);  xx.append(x[i]); yy.append(y[i])
i = 2  ;   annotation(note[i], x[i], y[i], x[i]+0.25, y[i]-0.11);  xx.append(x[i]); yy.append(y[i])
i = 3  ;   annotation(note[i], x[i], y[i], x[i]+0.30, y[i]+0.08);  xx.append(x[i]); yy.append(y[i])
i = 4  ;   annotation(note[i], x[i], y[i], x[i]+0.35, y[i]+0.035); xx.append(x[i]); yy.append(y[i])
i = 5  ;   annotation(note[i], x[i], y[i], x[i]+0.32, y[i]-0.075); xx.append(x[i]); yy.append(y[i])
i = 6  ;   annotation(note[i], x[i], y[i], x[i]+0.35, y[i]-0.03);  xx.append(x[i]); yy.append(y[i])
i = 7  ;   annotation(note[i], x[i], y[i], x[i]+0.30, y[i]-0.10);  xx.append(x[i]); yy.append(y[i])
i = 8  ;   annotation(note[i], x[i], y[i], x[i]+0.23, y[i]-0.13);  xx.append(x[i]); yy.append(y[i])
i = 9  ;   annotation(note[i], x[i], y[i], x[i]+0.35, y[i]+0.07);  xx.append(x[i]); yy.append(y[i])
i = 10 ;   annotation(note[i], x[i], y[i], x[i]-0.14, y[i]+0.07);  xx.append(x[i]); yy.append(y[i])
i = 11 ;   annotation(note[i], x[i], y[i], x[i]+0.35, y[i]-0.05);  xx.append(x[i]); yy.append(y[i])
i = 12 ;   annotation(note[i], x[i], y[i], x[i]-0.12, y[i]+0.05);  xx.append(x[i]); yy.append(y[i])
i = 13 ;   annotation(note[i], x[i], y[i], x[i]-0.30, y[i]);       xx.append(x[i]); yy.append(y[i])
i = 14 ;   annotation(note[i], x[i], y[i], x[i]-0.02, y[i]+0.12);  xx.append(x[i]); yy.append(y[i])
#i = 15 ;   annotation(str(i) + note[i], x[i], y[i], x[i]-0.20, y[i]-0.05)
#i = 16 ;   annotation(str(i) + note[i], x[i], y[i], x[i], y[i])
i = 17 ;   annotation(note[i], x[i], y[i], x[i]-0.30, y[i]-0.05);  xx.append(x[i]); yy.append(y[i])
i = 18 ;   annotation(note[i], x[i], y[i], x[i]-0.20, y[i]+0.00);  xx.append(x[i]); yy.append(y[i])
i = 19 ;   annotation(note[i], x[i], y[i], x[i]+0.18, y[i]+0.00);  xx.append(x[i]); yy.append(y[i])
#i = 20 ;   annotation(str(i) + note[i], x[i], y[i], x[i], y[i])
#i = 21 ;   annotation(str(i) + note[i], x[i], y[i], x[i], y[i])
#i = 22 ;   annotation(str(i) + note[i], x[i], y[i], x[i], y[i])
#i = 23 ;   annotation(str(i) + note[i], x[i], y[i], x[i], y[i]);
i = 24 ;   annotation(note[i], x[i], y[i], x[i]-0.10, y[i]+0.08);  xx.append(x[i]); yy.append(y[i])

plt.plot([0, 3], [0, 3], "k-", lw=1.5)
plt.plot(xx, yy, "bo")
plt.plot(x, y2, "ro")

Display("$\\gamma$ ($eV \cdot atom^{\minus 1}$)", "Predicted $\\gamma$ ($eV \cdot atom^{\minus 1}$)", [min(x)-min(x)*0.1, max(x)*1.02],
        [min(x)-min(x)*0.1, max(x)*1.02], "")

