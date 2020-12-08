



import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


site = "dist_" + str(sys.argv[2])                                         # site description, e.g. Mg
name = sys.argv[1][:-4]
data = np.loadtxt(sys.argv[1], comments="#")                             # import and reads data
max_columns = len(data[1])
column_labels = open(sys.argv[1], "r").readlines()[max_columns + 7]         # labels are the firts raw of data
label = column_labels.split()
label.pop(0)

for i in range(len(label)):
    if label[i] == site:
        x = data[:, i]
    elif label[i] == "Eadh":
        y = data[:, i]

def morse(x, d_eq, a, r_eq):
    return d_eq * (np.exp(-2 * a * (x - r_eq)) - 2 * np.exp(-a * (x - r_eq)))     # MORSE potential

for i in range(len(y)):
    if y[i] == min(y):
        r0_eq = x[i]
        d0_eq = np.abs(y[i])

popt, pcov = curve_fit(morse, x, y, bounds=([d0_eq*0.8, 0, r0_eq*0.8], [d0_eq*1.2, 50, r0_eq*1.2]))
d_eq, a, r_eq = popt
r2 = 1-np.sqrt(sum([(y[i] - morse(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
print(popt, r2)
out = open("MP_trendline", "w+")
out.write("d_eq= {} a= {} r_eq= {} R^2= {}" .format(d_eq, a, r_eq, r2))
out.close()
label_trend = "$E_{Adh} = %.3f \cdot e^{\minus %.3f (r - %.3f)} \minus 2 \cdot " \
              "e^{\minus %.3f (r \minus %.3f)}$ ; $R^{2}$= %.2f" %(d_eq, 2*a, r_eq, a, r_eq, r2)

if sys.argv[2] == "O":
    x_min = 1.5
    x_max = 5
    y_min = -3
    y_max = 0.1
elif sys.argv[2] == "Mg":
    x_min = 2
    x_max = 5
    y_min = -3
    y_max = 0.1

plt.plot(x, y, "ko") #, label=str(sys.argv[2]))
plt.plot([0, x_max], [0, 0], "k--") #, label=str(sys.argv[2]))
plt.plot(np.linspace(x_min, x_max, 101), morse(np.linspace(x_min, x_max, 101), *popt),
         "b:", label=label_trend)

plt.xlabel("$dist_{" + str(sys.argv[2]) + "}$" +" ($\\AA$)", fontsize=14)
plt.ylabel("$E_{Adh}$ (eV)", fontsize=14)
#plt.xticks(np.arange(0,Xmax, 0.25))
plt.xlim([x_min, x_max])#[min(x)*0.9, max(x)*1.2])
plt.ylim([y_min, y_max]) #[min(y)*1.2, 0.1])
#plt.tick_params(axis='both', labelrotation=0, labelsize=16)               # custimise tick labels
#plt.grid(True)
plt.legend(loc='best')
plt.title("$TM_{" + str(int(data[0,0])) + "}$ $dist_{" + str(sys.argv[2]) + "}$")
plt.savefig("MP" + str(int(data[0,0])) + "_" + site +".png", figsize=(11.69, 16.53), clear=True, dpi=300, orientation='landscape', transparent=True)
plt.show()


