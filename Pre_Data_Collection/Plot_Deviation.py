'''

	USAGE: ~.py input.dat
   	input: x,y file with comments for labels

'''

import sys
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt


icolour = ["b", "r", "k", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] ## n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))]


def SaveFig():
	answer = str(input("Would you like to save the figure (y/n)?\n"))
	if answer == "y":
		figure_out_name = str(input("What would it be the figure name (a word & no format)?\n"))
		plt.savefig(figure_out_name + ".svg", figsize=(12, 10), clear=True,
										bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)

ifile = open(sys.argv[1]).readlines()
x = {}
y = {}
for i in range(len(ifile)):
    line = ifile[i].split()
    if len(line) > 0:
        if line[1] == "system:":
            label = str(line[2] + "$^{\cdot c =" + line[4] + "} \cdot \\tau \leq$ " + line[6] + "eV")
        elif line[1].startswith("$E"):
            x_label = str(line[1] + " " + line[2] + " " + line[3] + " " + line[4] + " " + line[5])
            y_label = str(line[6] + " " + line[7] + " " + line[8] + " " + line[9] + " " + line[10] + " " + line[11])
        else:
            if label in x:
                x[label].append(float(line[0]))
                y[label].append(float(line[1]))
            else:
                x[label] = [float(line[0])]
                y[label] = [float(line[1])]

for n, label in enumerate(x):
    n_marker = n
    n_colour = n
    if n >= 2*len(icolour):
        n_colour = n - 2*len(icolour)
    elif n >= len(icolour):
        n_colour = n - len(icolour)
    if n >= len(imarker):
        n_marker = n - len(imarker)
    plt.plot(x[label], y[label], linestyle="None", marker=imarker[n], color=icolour[n], markersize=3, label=label)

x_lim = [min([min([i for i in x[j]]) for j in x])*1.1, max([max([i for i in x[j]]) for j in x])*1.1]
y_lim = [min([min([i for i in y[j]]) for j in y])*1.1, max([max([i for i in y[j]]) for j in y])*1.1]

plt.plot(x_lim, x_lim, "k-", lw=1.5)

plt.xlabel(x_label, fontsize=14)
plt.ylabel(y_label, fontsize=14)
plt.tick_params(axis='both', labelrotation=0, labelsize=12)               # custimise tick labels
plt.xlim(x_lim)
plt.ylim(y_lim)
plt.legend(loc='best')
plt.tight_layout()
#plt.ion()
plt.show()
SaveFig()
#plt.clf()
