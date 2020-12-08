



import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


data = np.loadtxt(sys.argv[1], comments="#")                    # import and reads data
depth_data= np.loadtxt(sys.argv[2], comments="#")                    # import and reads data

x = data[:, 1] * data[:, 4]
y = depth_data[:, 0]
xlabel = "$c_{i} \cdot i_{cc}$"
ylabel = "Depth /eV"

def pol(x, a, b, c):                                            # LOGARITHMIC
    return a + b * np.log(x + c)

popt, pcov = curve_fit(pol, x, y, bounds=([-50, -50, 0], [50, 50, 50]))
a, b, c = popt
r2 = 1-np.sqrt(sum([(y[i] - pol(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
trend_label = "$Depth = %.3f \plus %.3f \cdot ln(c_{i} \cdot i_{cc} \plus %.3f)$; $R^{2} = %.2f $" %(a, b, c, r2)
print(popt, r2)

plt.plot(np.linspace(min(x)*0.9, max(x)*1.1, 101), pol(np.linspace(min(x)*0.9, max(x)*1.1, 101), *popt),
         "b:", lw=1.5, label=trend_label)
plt.plot(x, y, "ko", markersize=3) #, label="column_" + str(y_column))
plt.xlabel(xlabel, fontsize=16)
plt.ylabel(ylabel, fontsize=16)
plt.tick_params(axis='both', labelrotation=0, labelsize=14)    		# custimise tick labels
plt.legend(loc='best')
plt.tight_layout()
#plt.grid(True)
#plt.savefig("Area.png", figsize=(11.69,16.53), clear=True, dpi=300, orientation='portrait',transparent=True)
plt.show()


