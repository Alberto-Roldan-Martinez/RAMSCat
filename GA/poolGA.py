'''
Birmingham Parallel Genetic Algorithm

A pool genetic algorithm for the
structural characterisation of 
nanoalloys.

Please cite - 
A. Shayeghi et al, PCCP, 2015, 17, 2104-2112

Authors -
The Johnston Group

13/2/15


--- Main PoolGA Class ---  

'''

import os
from random import randrange, uniform

from MinimisePool import minPool
from MinimiseOff import minOff
from MinimiseMut import minMut
from checkPool import checkPool as checkPool

import sys

class poolGA:

	def __init__(self,natoms,r_ij
				,eleNums,eleNames
				,eleMasses,mutate
				,nPool,cross,mutType
				,subString
				,boxAdd
				,surface
				,surfGA):
		
		self.nPool = nPool
		self.r_ij = r_ij
		self.mutate = mutate
		self.natoms = natoms
		self.eleNums = eleNums
		self.eleNames = eleNames
		self.eleMasses = eleMasses
		self.cross = cross
		self.mutType = mutType
		self.subString = subString
		self.boxAdd = boxAdd

		'''
		Surface Object.
		'''

		self.surface = surface
		self.surfGA = surfGA

		''' --- ''' 

		self.stride = natoms + 2

		self.runJobs()

	def runJobs(self):
		print("   poolGA")

		notFinished = self.checkFinished()

		while notFinished:

			pool = minPool(self.natoms,self.r_ij
					,self.eleNums,self.eleNames
					,self.eleMasses,self.nPool
					,self.stride,self.subString
					,self.boxAdd
					,self.surface,self.surfGA)

			notFinished = self.checkFinished()

		for i in range(self.nPool, self.nPool + 10*(self.natoms - 1)):			# Alberto 31/05/2022: increase from 1000 to 500*(atoms-1) >>>03/10/2022 changed to 100*(natoms)
			print("   poolGA", "i=", i, "i_max=", 10*sum(self.natoms-1))

			check = checkPool()
			converged = check.Convergence()

			if self.checkRunning():

				''' 
				If no new strucs are available
				and the some are still running
				allow the pool to grow.
				'''

				off = minMut(self.natoms,self.r_ij
					,"random",self.eleNums
					,self.eleNames,self.eleMasses
					,self.nPool,self.stride
					,self.subString,self.boxAdd
					,self.surface,self.surfGA)

			else:

				'''
				Else if pool is fully 
				minimimised begin mutation 
				and crossover.
				''' 

				self.decide()

	def decide(self):

		mutateRate = self.mutate * self.nPool

		choice = uniform(0,self.nPool)

		if choice < mutateRate:

			off = minMut(self.natoms,self.r_ij
				,self.mutType,self.eleNums
				,self.eleNames,self.eleMasses
				,self.nPool,self.stride
				,self.subString,self.boxAdd
				,self.surface,self.surfGA)
		else:

			off = minOff(self.natoms,self.eleNums
				,self.eleNames,self.eleMasses
				,self.nPool,self.cross,self.stride
				,self.subString,self.boxAdd
				,self.surface,self.surfGA)

	def checkRunning(self):

		'''
		Check if pool
		is running.
		'''

		with open("pool.dat","r") as pool:
			for line in pool:
				if "Running" in line:
					return True
				elif "Restart" in line:
					return True

		return False

	def checkFinished(self):

		'''
		Checks if all initial
		geometries relaxed.
		'''

		with open("pool.dat","r") as pool:
			for line in pool:
				if "NotMinimised" in line:
					return True
				elif "Restart" in line:
					return True

		return False
