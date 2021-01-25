



import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


data_coord = np.loadtxt(sys.argv[1], comments="#")                    # import and reads data
data = open(sys.argv[1]).readlines()

x0 = data_coord[:, 0]
y0 = [i*1E20 for i in data_coord[:, 1]]         # 1= Areas in m^2
#y0 = data_coord[:, 2]                           # 2= SurfEnergy


def SaveFig():
        answer = str(input("Would you like to save the figure (y/n)?\n"))
        if answer == "y":
                figure_out_name = str(input("What would it be the figure name (a word & no format)?\n"))
                plt.savefig(figure_out_name + ".svg", figsize=(12, 10), clear=True,
                                        bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)

def annotation(note, x0, y0, x, y):
    return  plt.annotate(note, xy=(x0,y0), xycoords="data", size=14,
                               xytext=(x,y), textcoords="data",
                               arrowprops=dict(arrowstyle="<-", color="b", fc="0.6"),
                               horizontalalignment="right", verticalalignment="top")
def notes(i):
    line = DataSet[i]
    if line[4] != str(0):
        note = str("$\it p$("+line[2]+")$\cdot$("+line[3]+")"+line[4])
    else:
        note = str("$\it p$("+line[2]+")$\cdot$("+line[3]+")")
    return note

def trendline(x, a, b, c, d):
#    return a*x + b                      # linear
#    return a*np.cos(b+x*c)
#    return (a / (np.sqrt(2*np.pi) * c )) * np.exp(-(x-b)**2 / (2*c**2))    # GAUSSIAN
    return a + (2 * -c * d/np.pi) / (4 * (x - b)**2 - c**2)     # LORENTZIAN


# Find the Trendline
popt, pcov = curve_fit(trendline, x, y) #, bounds=([d0_eq*0.8, 0, r0_eq*0.8], [d0_eq*1.2, 50, r0_eq*1.2]))
r2 = 1-np.sqrt(sum([(y[i] - trendline(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
if len(popt) == 2:
    a, b = popt
    trend_label = "$Area = %.3f %s + %.3f ; R^{2}$= %.2f$" %(a, "coord", b, r2)
elif len(popt) == 3:
    a, b, c = popt
    trend_label = "$Area = %.3f cos(%.3f \plus %.3f %s ; R^{2}$= %.2f$" %(a, b, c, "coord", r2)     # cos
elif len(popt) == 4:
    a, b, c, d = popt
    trend_label = "$Area = %.3f - \frac{%.3f}{4 * (%s - %.3f) - %.3f} ; R^{2}$= %.2f$" %(a, 2*c*d/np.pi, "coord", b**2, c**2, r2)

r2 = 1-np.sqrt(sum([(y[i] - trendline(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
print(popt, r2)

figure1 = plt.figure(figsize=(12, 10), clear=True)
plt.plot(x0, y0, "o")
plt.plot(np.linspace(0, 12, 150), trendline(np.linspace(0, 12, 150), *popt), "b--", label=trend_label)
plt.xlabel("Coordination", fontsize=14)
plt.ylabel("Area $(\\AA^{2} \cdot atom^{\minus 1})$", fontsize=14)
#       plt.xticks(np.arange(0,Xmax, 0.25))
plt.xlim([0, 12])
plt.ylim([0, max(y)*1.1])
plt.tick_params(axis='both', labelrotation=0, labelsize=12)               # custimise tick labels
#plt.title("$TM_{" + str(int(data[0,0])) + "}$ $dist_{" + str(sys.argv[2]) + "}$")
plt.legend(loc='best')
plt.grid(True)
plt.tight_layout()
plt.ion()
plt.show()
SaveFig()








# What to plot?
x = [float(DataSet[i][1]) for i in range(len(DataSet)) if float(DataSet[i][0]) > 0]
Y = [float(DataSet[i][1]) for i in range(len(DataSet)) if float(DataSet[i][0]) > 0]

plt.plot(np.linspace(0,max(X)+10,2),np.linspace(0,max(X)+10,2), "k-", lw=1.5)
plt.plot(X, Y,"rs", markersize=7)
#plt.plot(X2, Y2,"rs", fillstyle = "none", markersize=7)


plt.xlabel("Area /$\\AA ^{2}$", fontsize=16)
plt.ylabel("Predicted Area /$\\AA ^{2}$", fontsize=16)
plt.tick_params(axis='both',labelrotation=0,labelsize=14)           # custimise tick labels

i = 0  ;  x = X[i]+50  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 1  ;  x = X[i]+50  ;  y = Y[i]+5   ; annotation(notes(i),X[i],Y[i],x,y)
i = 2  ;  x = X[i]+75  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 3  ;  x = X[i]+75  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 4  ;  x = X[i]+85  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 5  ;  x = X[i]+60  ;  y = Y[i]+5   ; annotation(notes(i),X[i],Y[i],x,y)
i = 6  ;  x = X[i]-20  ;  y = Y[i]+10  ; annotation(notes(i),X[i],Y[i],x,y)
i = 7  ;  x = X[i]-15  ;  y = Y[i]+15  ; annotation(notes(i),X[i],Y[i],x,y)
i = 8  ;  x = X[i]+70  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 9  ;  x = X[i]-15  ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 10 ;  x = X[i]+140 ;  y = Y[i]+10  ; annotation(notes(i),X[i],Y[i],x,y)
i = 11 ;  x = X[i]-35  ;  y = Y[i]+20  ; annotation(notes(i),X[i],Y[i],x,y)
i = 12 ;  x = X[i]+60  ;  y = Y[i]+5   ; annotation(notes(i),X[i],Y[i],x,y)
i = 13 ;  x = X[i]+80  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 14 ;  x = X[i]+120 ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 15 ;  x = X[i]+60  ;  y = Y[i]+8   ; annotation(notes(i),X[i],Y[i],x,y)
i = 16 ;  x = X[i]+100 ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 17 ;  x = X[i]+50  ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 18 ;  x = X[i]-5   ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 19 ;  x = X[i]-30  ;  y = Y[i]+20  ; annotation(notes(i),X[i],Y[i],x,y)
i = 20 ;  x = X[i]-10  ;  y = Y[i]+2   ; annotation(notes(i),X[i],Y[i],x,y)
i = 21 ;  x = X[i]-30  ;  y = Y[i]+5   ; annotation(notes(i),X[i],Y[i],x,y)
i = 22 ;  x = X[i]+20  ;  y = Y[i]+30  ; annotation(notes(i),X[i],Y[i],x,y)
i = 23 ;  x = X[i]-25  ;  y = Y[i]+30  ; annotation(notes(i),X[i],Y[i],x,y)
i = 24 ;  x = X[i]-20  ;  y = Y[i]+33  ; annotation(notes(i),X[i],Y[i],x,y)
i = 25 ;  x = X[i]-15  ;  y = Y[i]+40  ; annotation(notes(i),X[i],Y[i],x,y)
i = 26 ;  x = X[i]-18  ;  y = Y[i]+23  ; annotation(notes(i),X[i],Y[i],x,y)
i = 27 ;  x = X[i]+80  ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 28 ;  x = X[i]+50  ;  y = Y[i]+5   ; annotation(notes(i),X[i],Y[i],x,y)
i = 29 ;  x = X[i]-20  ;  y = Y[i]+22  ; annotation(notes(i),X[i],Y[i],x,y)
i = 30 ;  x = X[i]-20  ;  y = Y[i]     ; annotation(notes(i),X[i],Y[i],x,y)
i = 31 ;  x = X[i]+55  ;  y = Y[i]+3   ; annotation(notes(i),X[i],Y[i],x,y)
i = 32 ;  x = X[i]-20  ;  y = Y[i]+20  ; annotation(notes(i),X[i],Y[i],x,y)


plt.grid(True)
plt.tight_layout()
plt.ion()
plt.show()
SaveFig()
