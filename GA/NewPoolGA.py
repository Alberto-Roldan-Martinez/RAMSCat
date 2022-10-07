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

--- BPGA Input file ---

'''

import os
from random import uniform, randint

from GA.MinimiseRan import minRan 
from GA.MinimiseOff import minOff
from GA.MinimiseMut import minMut
import GA.Database as db						# Alberto 3/10/2022 defining the number of calculations
#from GA.checkPool import checkPool as checkPool

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

		self.run()

	def run(self):

		'''
		Minimise random geometries and
		add them to the pool.dat file.
		'''
		while self.getPoolSize() < self.nPool:
			print("   Adding to pool ", db.findLastDir() + 1, " - NewPoolGA")
			pool = minRan(self.natoms,self.r_ij
						  ,self.eleNums,self.eleNames
						  ,self.eleMasses,self.nPool
						  ,self.stride,self.subString
						  ,self.boxAdd
						  ,self.surface,self.surfGA)
		'''
		Once nPool structure are in
		pool.dat begin crossover
		and mutation. 
		'''
		'''
		Set the mutation rate.
		'''
		mutateRate = self.mutate * self.nPool

#		for i in range(sum(self.eleNums)*100):				# Alberto 01/08/2022: Changed 1000 by sum(eleNums)*100 --> eleNums is list >> 03/10/2022 change "for" with "while"
		while db.findLastDir() < sum(self.eleNums)*5:	#####################
			choice = uniform(0, self.nPool)
			if choice < mutateRate:
				mutantChoice = randint(0, len(self.mutType) - 1)		# Alberto 07/10/2022 Various mutant types added at once, randomly chosen at every mutant step
				print("   Mutant:", self.mutType[mutantChoice], "structure", db.findLastDir() + 1, " - NewPoolGA")
				off = minMut(self.natoms,self.r_ij
					,self.mutType[mutantChoice], self.eleNums
					,self.eleNames,self.eleMasses
					,self.nPool,self.stride
					,self.subString,self.boxAdd
					,self.surface,self.surfGA)
			else:
				print("   Offspring", db.findLastDir() + 1, " - NewPoolGA")
				off = minOff(self.natoms,self.eleNums
					,self.eleNames,self.eleMasses
					,self.nPool,self.cross,self.stride
					,self.subString,self.boxAdd
					,self.surface,self.surfGA)

	def getPoolSize(self):

		if os.path.exists("pool.dat") == False:
			return 0
		else: 
			with open("pool.dat","r") as pool:
				poolList = pool.readlines()
				poolSize = len(poolList) / (self.natoms + 2)
				return poolSize 

