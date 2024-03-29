'''
Birmingham Parallel Genetic Algorithm

A pool genetic algorithm for the
structural characterisation of 
nanoalloys.

Please cite - 
A. Shayeghi et al, PCCP, 2015, 17, 2104-2112

Authors -
The Johnston Group

16/10/14

--- checkPool Class ---

'''

class checkPool:

    def __init__(self):

        self.energies = []
        self.getEnergies()

    def getEnergies(self):

        '''
        Get final energies
        from pool.dat after
        convergence.
        '''

        with open("pool.dat","r") as pool:
            for line in pool:
                if "Energy" in line:
                    energyList = line.split()
                    self.energies.append(float(energyList[2]))
                elif "Running" in line:
                    self.energies.append(-1000000)
					
    def checkEnergy(self, newEnergy, sphericity):

        HighestEnergy = max(self.energies)
        self.lowestIndex = self.energies.index(HighestEnergy)

        if newEnergy < HighestEnergy and sphericity <= 1.:              # Alberto 05/08/2022
            self.Index = self.energies.index(HighestEnergy)
            return True
        else:
            return False

    def Convergence(self):

        poolMin = min(self.energies)
        poolMax = max(self.energies)

        if poolMin - poolMax > -0.01:
            print ("Calculation Converged")
            return True
        else:
            return False
