'''
Birmingham Parallel Genetic Algorithm

A pool genetic algorithm for the
structural characterisation of 
nanoalloys.

Please cite - 
A. Shayeghi et al, PCCP, 2015, 17, 2104-2112

Authors -
The Johnston Group

20/3/15

--- Offspring Minimiser Class ---

'''

import os

import GA.Database as db

from GA.DFT_input import vasp_input as DFTin
from GA.DFT_output import vasp_output as DFTout

from GA.Crossover import crossover as cross 
from GA.checkPool import checkPool as checkPool
from GA.CoM import CoM 
from GA.Explode import checkClus
from GA.fixOverlap import fixOverlap

from GA.SurfOpt import SurfOpt 
from GA.surfacePOSCAR import surfacePOSCAR 

import sys

class minOff: 

	def __init__(self,natoms,eleNums
				,eleNames,eleMasses
				,nPool,cross,stride
				,subString,boxAdd
				,surface,surfGA):
		
		self.natoms = natoms
		self.eleNames = eleNames
		self.eleMasses = eleMasses
		self.eleNums = eleNums
		self.nPool = nPool
		self.cross = cross
		self.stride = stride
		self.subString = subString
		self.boxAdd = boxAdd

		'''
		Surface Object.
		'''

		self.surface = surface
		self.surfGA = surfGA

		''' --- ''' 
		
		self.runCalc()

	def runCalc(self):

		'''
		Start calculation
		making new 
		directory.
		'''

		db.lock()

		self.calcNum = db.findLastDir() + 1

		while os.path.exists(str(self.calcNum)): 
			self.calcNum = db.findLastDir() + 1 

		os.system("mkdir " + str(self.calcNum))

		self.produceOffspring()

		db.unlock()

		self.minimise()

	def produceOffspring(self):

		'''
		Produces XYZ after 
		crossover. 
		'''

		newClus = cross(self.cross,self.nPool,self.stride
						,self.eleNums,self.eleNames,self.natoms)

		self.offspring = newClus.mate()
		self.offspring = fixOverlap(self.offspring)

		if self.surfGA:

			SurfaceStruc = SurfOpt(self.offspring,self.surface,self.eleNames,self.eleMasses)

			SurfClus = SurfaceStruc.placeClus()

			self.vaspIN = surfacePOSCAR(self.calcNum,self.eleNames,SurfClus,self.surface)

		else:

			self.vaspIN = DFTin(self.calcNum,self.offspring,self.eleNames
								,self.eleMasses,self.eleNums,self.boxAdd)

	def restart(self):

		'''
		Restart Calculation
		without making 
		new directory.
		'''

		db.lock()

		self.produceOffspring()

		db.unlock()

		self.minimise()

	def minimise(self):

		'''
		Start 
		DFT calculation.
		'''

		if self.doDFT() == 0:

			output = DFTout(self.calcNum,self.natoms)

			if output.checkError():
				self.restart()
			else:
				self.finalEnergy, self.sphericity = output.getEnergy()
				self.finalCoords = output.getCoords()

				check = checkClus(self.natoms,self.finalCoords)

				if check.exploded() == False:
					self.updatePool()
				else:
					self.restart()

	def doDFT(self):

		'''
		Change directory and 
		submit calculation.
		'''

		base = os.environ["PWD"]
		os.chdir(base+"/"+str(self.calcNum))

		exitcode = os.system(self.subString)

		with open(base+"/exitcodes.dat","a") as exit:
			exit.write(str(self.calcNum))
			exit.write(" Exitcode = "+str(exitcode)+" Offspring\n")
			
		os.chdir(base)

		return exitcode

	def updatePool(self):

		AcceptReject = checkPool()
		Accept = AcceptReject.checkEnergy(float(self.finalEnergy), float(self.sphericity))

		if Accept:
			index = AcceptReject.lowestIndex
			index = (index*self.stride)+1

			db.updatePool("Finish"
					,self.calcNum
					,index,self.eleNums
					,self.eleNames,self.eleMasses
					,self.finalEnergy, self.sphericity, self.finalCoords
					,self.stride,self.vaspIN.box)
			print("   Offspring:", db.findLastDir() + 1)