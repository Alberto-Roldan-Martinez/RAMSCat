



import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

data_coord = np.loadtxt(sys.argv[1], comments="#")                    # import and reads data
data = open(sys.argv[2]).readlines()
data = [data[i].split() for i in range(len(data))]

x0 = data_coord[:, 0]
y0 = data_coord[:, 2]                           # 2= SurfEnergy


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
                               horizontalalignment="right", verticalalignment="center")

def trendline(x, a, b, c, d):
#    return a*x + b                      # linear
#    return a*np.cos(b+x*c)
#    return (a / (np.sqrt(2*np.pi) * c )) * np.exp(-(x-b)**2 / (2*c**2))    # GAUSSIAN
#    return a + (2 * -c * d/np.pi) / (4 * (x - b)**2 - c**2)     # LORENTZIAN
    return a - (b + x**c)/d**2      # LORENTZIAN


# Find the Trendline
popt, pcov = curve_fit(trendline, x0, y0)   #, bounds=(-50, [150, 150, 150, 100]))
r2 = 1-np.sqrt(sum([(y0[i] - trendline(x0[i], *popt))**2 for i in range(len(x0))])/sum(i*i for i in y0))
a, b, c, d = popt
#trend_label = str("$\\gamma$ = %.3f %s + %.3f ; R^{2}$= %.2f$" %(a, "coord", b, r2))
#trend_label = "$\\gamma$ = %.3f cos(%.3f \plus %.3f %s ; R^{2}$= %.2f$" %(a, b, c, "coord", r2)     # cos
trend_label = "Lorentzian: a, b, c, d =" + str(round(a, 3)) + ", " +\
                  str(round(b, 3)) + ", " + str(round(c, 3)) + ", " + str(round(d, 3))

plt.plot(x0, y0, "bo")
plt.plot(np.linspace(0, 12, 150), trendline(np.linspace(0, 12, 150), *popt), "b--", label="$R^{2}$= " + str(round(r2, 2)))

plt.xlabel("Coordination", fontsize=14)
plt.ylabel("$\\gamma$ ($J \cdot m^{\minus 2}$)", fontsize=14)
plt.tick_params(axis='both', labelrotation=0, labelsize=12)               # custimise tick labels
plt.xlim([0, 12.15])
plt.ylim([min(y0)-0.1, max(y0)*1.15])
plt.title(trend_label)
plt.legend(loc='best')
#plt.grid(True)
#plt.tight_layout()
#plt.ion()
plt.show()
#SaveFig()
plt.clf()


# Plot the real vs. estimated
# What to plot?
coord = [float(data[i][4]) for i in range(len(data))]               # coordination to get predicted area
natoms = [float(data[i][3]) for i in range(len(data))]              # number of atoms with such coordination

x = [float(data[i][2]) for i in range(len(data))]              # real Area
y = [trendline(coord[i], *popt) for i in range(len(coord))]

print(x, "\n", y)

plt.plot([min(x), max(x)], [min(x), max(x)], "k-", lw=1.5)
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

i = 0  ;   annotation(str(i) + note[i], x[i], y[i], x[i], y[i])
#i = 1  ;   annotation(str(i) + note[i], x[i], y[i], x[i]+300, y[i]-55)
#i = 2  ;   annotation(str(i) + note[i], x[i], y[i], x[i]+300, y[i]+0)
#i = 3  ;   annotation(str(i) + note[i], x[i], y[i], x[i]+300, y[i]-75)
#i = 4  ;   annotation(str(i) + note[i], x[i], y[i], x[i]+310, y[i]-35)
#i = 5  ;   annotation(str(i) + note[i], x[i], y[i], x[i]+300, y[i]+40)
#i = 6  ;   annotation(str(i) + note[i], x[i], y[i], x[i]-80, y[i]+5)
#i = 7  ;   annotation(str(i) + note[i], x[i], y[i], x[i]+180, y[i]-0)
#i = 8  ;   annotation(str(i) + note[i], x[i], y[i], x[i]-100, y[i]+45)
#i = 9  ;   annotation(str(i) + note[i], x[i], y[i], x[i]-80, y[i]+15)
#i = 10 ;   annotation(str(i) + note[i], x[i], y[i], x[i]-70, y[i]+25)
#i = 11 ;   annotation(str(i) + note[i], x[i], y[i], x[i]-70, y[i]+65)
#i = 12 ;   annotation(str(i) + note[i], x[i], y[i], x[i]-50, y[i]+40)
#i = 13 ;   annotation(str(i) + note[i], x[i], y[i], x[i]+130, y[i]-40)
#i = 14 ;   annotation(str(i) + note[i], x[i], y[i], x[i]+300, y[i]+17)
#i = 15 ;   annotation(str(i) + note[i], x[i], y[i], x[i]+15, y[i]+40)
#i = 16 ;   annotation(str(i) + note[i], x[i], y[i], x[i]+300, y[i]-20)
#i = 17 ;   annotation(str(i) + note[i], x[i], y[i], x[i]+180, y[i]+20)
#i = 18 ;   annotation(str(i) + note[i], x[i], y[i], x[i]-60, y[i]+80)
#i = 19 ;   annotation(str(i) + note[i], x[i], y[i], x[i]+150, y[i]+25)
#i = 20 ;   annotation(str(i) + note[i], x[i], y[i], x[i]+250, y[i]+18)
#i = 21 ;   annotation(str(i) + note[i], x[i], y[i], x[i]+250, y[i]+40)
#i = 22 ;   annotation(str(i) + note[i], x[i], y[i], x[i]-200, y[i]-20)
#i = 23 ;   annotation(str(i) + note[i], x[i], y[i], x[i]-200, y[i]-0)
#i = 24 ;   annotation(str(i) + note[i], x[i], y[i], x[i]+10, y[i]+50)


plt.xlabel("$\\gamma$ ($J \cdot m^{\minus 2}$)", fontsize=14)
plt.ylabel("Predicted $\\gamma$ ($J \cdot m^{\minus 2}$)", fontsize=14)
plt.tick_params(axis='both', labelrotation=0, labelsize=12)           # custimise tick labels
plt.xlim([1.20, 1.70])
plt.ylim([1.20, 1.70])
#plt.title(trend_label)
#plt.legend(loc='best')
#plt.grid(True)
plt.tight_layout()
#plt.ion()
plt.show()
#SaveFig()







