



import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


DataSet = np.loadtxt(sys.argv[1], comments="#")                    # import and reads data
columnLabels = open(sys.argv[1],"r").readlines()[0]                  # labels are the firts raw of data
label = columnLabels.split()
label.pop(0)

DataSet_2 = np.loadtxt(sys.argv[2], comments="#")                    # import and reads data

# What to plot?
X = DataSet[:,0]
Y = DataSet[:,1]
Y2 = DataSet_2[:,1]


def function_A(x,a,b):                                            # x= coordination
    return a*x + b
def function_B(x,a,b,c,d):                                            # x= coordination
#    return a*np.cos(b+x*c)
#    return (a / (np.sqrt(2*np.pi) * c )) * np.exp(-(x-b)**2 / (2*c**2))    # GAUSSIAN
    return a + ( 2 * -c * d/np.pi) / (4 * ( x - b )**2 - c**2)     # LORENTZIAN


popt_A,pcov_A=curve_fit(function_A,X,Y)
perr=np.sqrt(np.diag(pcov_A))
Serror=sum(i*i for i in (Y - function_A(X, *popt_A)))
Stotal=sum(i*i for i in Y)
Root_squared_A=(1-np.sqrt(Serror/Stotal))*100

popt_B,pcov_B=curve_fit(function_B,X,Y2 , bounds=([0,0,0,0],[50,50,50,10000]))
perr=np.sqrt(np.diag(pcov_B))
Serror=sum(i*i for i in (Y2 - function_B(X, *popt_B)))
Stotal=sum(i*i for i in Y2)
Root_squared_B=(1-np.sqrt(Serror/Stotal))*100




fig, ax1 = plt.subplots()
ax1.set_xlim((0,12.25))
ax1.set_xlabel(label[0], fontsize=18)
ax1.set_ylabel("$\\gamma$ /$J \cdot m^{\minus 2}$", fontsize=18)
ax1.yaxis.label.set_color("b")
ax1.plot(X, Y ,"bo",markersize=7)
ax1.plot(np.linspace(0.5,12,101), function_A(np.linspace(0.5,12,101), *popt_A),"b:", label= "$\\gamma$ = %.3f * %s + %.3f ; $R^{2}$= %.2f %%" %(popt_A[0],label[0],popt_A[1],Root_squared_A))
ax1.tick_params(axis='both',labelrotation=0,labelsize=16)    # custimise tick labels
legend = ax1.legend(loc='best')
fig.tight_layout()
plt.grid(True)
plt.savefig("Energy_Trend.pdf", figsize=(11.69,16.53), clear=True, dpi=300, orientation='portrait',transparent=True)


fig, axB = plt.subplots()
axB.set_xlim((0,12.25))
axB.set_xlabel(label[0], fontsize=18)
#axB = ax1.twinx()
axB.set_ylabel("Area /$\\AA ^{2} \cdot atom^{\minus 1}$", fontsize=18)
axB.yaxis.label.set_color("r")
axB.plot(X, Y2 ,"rs",markersize=7)
axB.plot(np.linspace(0.5,12,101), function_B(np.linspace(0.5,12,101), *popt_B),"r:", label= "Area = %.3f - %.3f /\n\t(4 * (%s - %.3f)$^{2}$ - %.3f$^{2}$) ; $R^{2}$= %.2f %%" %(popt_B[0],2*popt_B[2]*popt_B[3]/np.pi,label[0],popt_B[1],popt_B[2],Root_squared_B))
axB.tick_params(axis='both',labelrotation=0,labelsize=16)    # custimise tick labels
legend = axB.legend(loc='best')
fig.tight_layout()
plt.grid(True)
plt.savefig("Area_Trend.pdf", figsize=(11.69,16.53), clear=True, dpi=300, orientation='portrait',transparent=True)

plt.show()


