'''
Birmingham Parallel Genetic Algorithm

A pool genetic algorithm for the
structural characterisation of 
nanoalloys.

Please cite - 
A. Shayeghi et al, PCCP, 2015, 17, 2104-2112

Authors -
The Johnston Group

19/1/15

--- Surface POSCAR Class ---

'''

import os, sys

class surfacePOSCAR():

    def __init__(self
                , calcNum
                , eleNames
                , slabClusXYZ
                , surface):

        self.calcNum = str(calcNum)
        self.eleNames = eleNames
        self.eleNamesSurf = []

#        self.addMgO()
        self.addsurface()    # ARM 14/11/2018

        self.xyz = slabClusXYZ
        self.x = surface.x
        self.y = surface.y
        self.z = surface.z
        self.vac = surface.vac
        self.lat = surface.lat
        self.box = 0.
        self.eleNums = []

        self.getEleNums()

        self.printPOSCAR()

#    def addMgO(self):
#
#        for element in self.eleNames:
#            self.eleNamesSurf.append(element)
#
#        self.eleNamesSurf.append("O")
#        self.eleNamesSurf.append("Mg")

    def addsurface(self):      # ARM 14/11/2018

        for element in self.eleNames:
            self.eleNamesSurf.append(element)

        self.eleNamesSurf.append("O")
        self.eleNamesSurf.append("Mg")
        self.eleNamesSurf.append("C")

    def getEleNums(self):

        '''
        Create a new, different eleNums 
        list containing cluster and slab.
        '''
        self.neweleNamesSurf = []      # ARM 14/11/2018
        for element in self.eleNamesSurf:
            eleCount = 0
            for i in self.xyz:
                ele, x, y, z = i
                if ele == element:
                    eleCount += 1

            if eleCount != 0:          # ARM 14/11/2018
                self.eleNums.append(eleCount)
                self.neweleNamesSurf.append(element)
        self.eleNamesSurf = self.neweleNamesSurf
 

    def printPOSCAR(self):

        for element in self.eleNamesSurf:         # ARM 15/11/2018
            if element == "Mg":
                xLat = str(self.x*(self.lat/2))
                yLat = str(self.y*(self.lat/2))
                zLat = str(self.z*(self.lat/2)+(self.vac * (self.lat/2)))
                xline = xLat + " 0.0 0.0"
                yline = "0.0 " + yLat + " 0.0"
                zline = "0.0 0.0 "+ zLat
            elif element == "C":
                xLat = str(self.x*self.lat)
                y1Lat = str(self.y*(-self.lat/2))
                y2Lat = str(self.y*0.866025405*self.lat)
                zLat = str(self.z *(self.lat)+(self.vac*(self.lat))) 
                xline = xLat + " 0.0 0.0"
                yline = y1Lat + " " + y2Lat + " 0.0"
                zline = "0.0 0.0 " + zLat

        with open(self.calcNum+"/POSCAR", "w") as poscar:

            poscar.write("Test\n")
            poscar.write("  1.0\n")

            poscar.write("    " + xline + "\n")
            poscar.write("    " + yline + "\n")
            poscar.write("    " + zline + "\n")

            ''' Write elements '''

            for element in self.eleNamesSurf:

                poscar.write("  "+element)

            poscar.write("\n")

            for eleNum in self.eleNums:
                poscar.write("  "+str(eleNum))

            '''
            Selective dynamics so that 
            the cluster relaxes and the 
            surface remains fixed.
            '''

            poscar.write("\nSelective Dynamics\n")
            
            '''
            Use cartesian coordinates.
            '''

            poscar.write("Cartesian\n")

            '''
            Write coordinates in order of 
            element names. If the coordinates 
            belong to a cluster all the geometry 
            to relax. 
            '''

            for element in self.eleNamesSurf:
                for i in self.xyz:
                    ele, x, y, z = i
                    if ele == element:
                        if ele == "Mg" or ele == "O" or ele == "C":
                            line = " "+str(x)+" "+str(y)+" "+str(z)+" F F F\n"
                        else:
                            line = " "+str(x)+" "+str(y)+" "+str(z)+" T T T\n"
                        poscar.write(line)
                        
#        os.system("cp {POTCAR,INCAR} "+self.calcNum)       # commented by Alberto 01/10/2021
