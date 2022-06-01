'''
Birmingham Parallel Genetic Algorithm

A pool genetic algorithm for the
structural characterisation of 
nanoalloys.

Please cite - 
A. Shayeghi et al, PCCP, 2015, 17, 2104-2112

Authors -
The Johnston Group

8/6/15

--- Graphene Generator Class ---

'''

import numpy as np



class GrapheneV1():

    def __init__(self, x, y, z, vac, clusHeight):

        self.xyz = []

        self.x = x
        self.y = y
        self.z = z

        self.vac = vac
        self.clusHeight = clusHeight

        self.lat = 2.4686

        self.genX()
        self.genY()
        self.genZ()
        self.genV()



    def getSurf(self): 

        return self.xyz

    def genX(self):

        for i in range(self.x):     #how many atoms to initialise
                x = i * self.lat
                y = 0.0
                z = 0.0
                line = ["C", x, y, z]
                self.xyz.append(line)
                
                x = i* self.lat
                y = 1.4247906635
                z = 0
                line = ["C", x, y, z]
                self.xyz.append(line)


    def genY(self):

        natoms = len(self.xyz)        #show how many atoms need to duplicate among y axis

        for i in range(1, self.y):     #how many times we need to duplicate
            for j in range(natoms):  #recognise different types of atoms
                i=float(i)
                ele, x, y, z = self.xyz[j]
                x = x - i/2*self.lat
                y = y + i*2.137870312
                z = 0

                newline = ["C", x, y, z]
                self.xyz.append(newline)


    def genZ(self):
        self.xyz


    def genV(self):
        x_diff = [np.abs(i[1]/self.x*self.lat/2 - 1) for i in self.xyz]
        y_diff = [np.abs(i[2]/self.y*self.lat/2 - 1) for i in self.xyz]
        xyz_v = []
        for i in range(len(self.xyz)):
            if x_diff[i] != min(x_diff) or y_diff[i] != min(y_diff):
                xyz_v.append(self.xyz[i])
        self.xyz = xyz_v


