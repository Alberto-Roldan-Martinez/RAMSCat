



import os, sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from Library import ecoh_bulk


name = sys.argv[1][:-4]
element = sys.argv[2]
dataset = np.loadtxt(sys.argv[1], comments="#")                             # import and reads data
max_columns = len(dataset[1])

icolour = ["b","r","g","c","m","k"]

bulk_coh, bulk_coord = ecoh_bulk(element)

x = dataset[:,1]
y = dataset[:,2]/bulk_coh
xlabel = "cc"
ylabel = "$\\frac{E_{Coh}}{E_{Coh}^{Bulk}}$"

def interpolation(x, a):
    bulk_coh, bulk_coord = ecoh_bulk(element)
    return np.log(a)/np.log(a/(a+bulk_coord)) - 1/np.log(a/(a+bulk_coord)) * np.log(a+x)

#def interpolation(x, a, b):
#    return a*x+b

#sigma = np.ones(len(x))
#sigma[[-1]] = 0.1               # gives more weight to the last point (bulk)
popt, pcov = curve_fit(interpolation, x, y)#, sigma=sigma)#, bounds=(0,5))  #bulk_coh, abs(bulk_coh)))
r2 = 1-np.sqrt(sum([(y[i] - interpolation(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))

#trend_label = "$ \\frac{E_{Coh}}{E_{Coh}^{Bulk}} = \\frac{%.3f}{log \\frac{%.3f}{%.3f}}" \
#              " \minus \\frac{1}{log \\frac{%.3f}{%.3f}} \cdot log(%.3f \plus cc);  R^2 = %.2f$" \
#              %(popt, popt, popt+bulk_coord, popt, popt+bulk_coord, popt, r2)
trend_label = "$ \\frac{E_{Coh}}{E_{Coh}^{Bulk}} = %.3f \minus (%.3f) \cdot log(%.3f \plus cc);  R^2 = %.2f$" \
              %(popt/np.log(popt/bulk_coord), 1/(np.log(popt/bulk_coord)), popt, r2)

#trend_label = "$ \\frac{E_{n,Coh}}{E_{Coh}^{Bulk}} = %.3f \cdot cc \plus %.3f$;  $R^2 = %.2f$" %(popt[0], popt[1], r2)


print("param =", popt)

figure = plt.figure(clear=True)
ax = figure.add_subplot()

ax.plot(x[:15], y[:15], "bs", fillstyle='none', markersize=3)
ax.plot(x[16:], y[16:], "ko",  markersize=3)
ax.plot(np.linspace(0, 12, 101), interpolation(np.linspace(0, 12, 101), *popt), "k:", label=trend_label)

ax.set_xlabel(xlabel, fontsize=16)
ax.set_ylabel(ylabel, fontsize=18)
ax.set_xlim(-0.15, 12.15)
ax.set_ylim(-0.05, 1.05)
#ax.xaxis.set_tick_params(labelsize=16)
#ax.yaxis.set_tick_params(labelsize=16)
#ax.xaxis.set_major_locator(LinearLocator(5))
#ax.yaxis.set_major_locator(LinearLocator(5))
ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
legend = ax.legend(loc='best')#, fontsize=16)
#plt.grid(True)
plt.subplots_adjust(top=0.98, bottom=0.12, left=0.15, right=0.98)
plt.savefig("Fitted_Ecoh_"+element+".png", dpi=300, orientation='landscape', transparent=True)
plt.show()


