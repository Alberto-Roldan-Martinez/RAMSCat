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

--- Database Functions ---

1. updatePool - Add final geometry to pool.dat
2. addEle - Add elements to final geometry from VASP.
3. findLastDir - Search for last calculation number. 
4. readPool - Read pool.dat into list. 
5. writePool - Write updated pool to pool.dat
6. lock - lock pool.dat
7. unlock - unload pool.dat

'''

import os, time, errno, sys

from GA.CoM import CoM 

global fd

def updatePool(upType
			,calcNum
			,index,eleNums
			,eleNames,eleMasses
			,finalEn, sphericity, finalCoords
			,stride,box):

	'''
	After completion take
	final energy and coords
	from OUTCAR and Update
	poolList
	'''

	lock()

	poolList = readPool()

	coordsEle=addEle(upType,finalCoords
					,eleNums,eleNames
					,eleMasses,box)

	poolList[index+1:index+stride-1]=coordsEle

	if "Finish" in upType:
		poolList[index]="Energy = "+str(finalEn)+" Dir = "+str(calcNum)+" Sphericity = " + str(sphericity) + "\n"
	elif "Restart" in upType:
		poolList[index]="Restart\n"

	writePool(poolList)

	unlock()

def addEle(upType,coords
		,eleNum,eleNames
		,eleMasses,box):

	'''
	Adds element types 
	to final coordinates
	from OUTCAR. Removes 
	from centre of box
	'''

	clus=[]
	eleList=[]
	coordsEle=[]

	'''
	For VASP output in
	centre of bounding box.
	'''

	if "Finish" in upType:
		coords = [float(i) - box/2 for i in coords]

	'''
	Add elements to list
	of coordinates. 
	'''

	for i in range(len(eleNames)):
		for j in range(eleNum[i]):
			eleList.append(eleNames[i])

	'''
	Convert to clus format
	for CoM conversion. 
	'''

	for i in range(0,len(coords),3):
		atom = [eleList[int(i/3)],coords[int(i)],coords[int(i+1)],coords[int(i+2)]]
		clus.append(atom)

	clus = CoM(clus,eleNames,eleMasses)

	'''
	Convert to list of 
	strings for pool update.
	'''

	for i in range(len(clus)):
		ele,x,y,z = clus[i]
		newLine = ele+" "+str(x)+" "+str(y)+" "+str(z)+"\n"
		clus[i] = newLine

	return clus

def findLastDir():

	'''
	Finds directory
	containing last
	calculation.
	'''

	calcList = []
	dirList = os.listdir(".")

	for i in dirList:
		try:
			calcList.append(int(i))
		except ValueError:
			continue

	calcList = sorted(calcList)

	'''
	If list is empty.
	'''

	if not calcList:
		lastCalc = 0
	else:
		lastCalc = calcList[len(calcList)-1]

	return lastCalc

def readPool():

	'''
	Reads pool at beginning/
	throughout calculation.
	'''

	with open("pool.dat","r") as pool:
		poolList = pool.readlines()

	return poolList

def writePool(poolList):

	'''
	Writes pool to
	file after any
	changes.
	'''

	with open("pool.dat","w") as pool:
		for line in poolList:
			pool.write(line)

def lock():

	'''
	An atomic operation on *most* OSes and filesystems
	Quote from man for open(2):
	
	"On NFS, O_EXCL is supported only when using NFSv3 or later on
	 kernel 2.6 or later.  In NFS environments where O_EXCL support
	 is not provided, programs that rely on it for performing
	 locking tasks will contain a race condition."
	'''

	global fd

	while True:

		try:
			fd = os.open('lock.db', os.O_CREAT | os.O_EXCL | os.O_RDWR)
			return

		except OSError as oserr:
			if oserr.errno == errno.EEXIST:
				time.sleep(0.5)	# Sleep for 500ms between polls
			else:
				raise

def unlock():

	global fd

	if (fd != None):
		os.close(fd)
		os.remove('lock.db')
		fd = None
	else:
		raise Exception("Attempted database unlock when not holding lock.")
		
