


import sys
import numpy as np

dataset = np.loadtxt(sys.argv[1], comments="#")
max_columns = len(dataset[1])
column_labels = open(sys.argv[1], "r").readlines()[max_columns + 7]
label = column_labels.split()
label.pop(0)

param = []

for i in range(len(label)-1):
    param.append(sum(dataset[:, i])/len(dataset[:, i]))

output = open("Param.dat", "w+")
for p in param:
    output.write("%.3f " % p)
output.write("\n")
output.close()