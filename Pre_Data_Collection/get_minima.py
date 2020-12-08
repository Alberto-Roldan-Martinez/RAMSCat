

import sys
import numpy as np


data = np.loadtxt(sys.argv[1], comments="#")                    # import and reads data

e_adh = data[:, 15]
done = 0

for i, line in enumerate(data):
    if done == 0:
        if line[15] == min(e_adh):
            output = open("min_Param.dat", "a+")
            for j in data[i, :]:
                output.write(" {:7.3f}"  .format(float(j)))
            output.write("\n")
            output.close()
            done = 1


