#!/bin/bash


if [ -f "Cluster/POSCAR" ];
	then rm Cluster/POSCAR Cluster/O* Cluster/CONTCAR
fi

sed -n "1,5 p" CONTCAR >> Cluster/POSCAR 
sed -n "6 p" CONTCAR | awk '{print $1}' >> Cluster/POSCAR  
f=$(sed -n "7 p" CONTCAR | awk '{print $1}') # >> Cluster/POSCAR   
echo $f >> Cluster/POSCAR
n=$(echo "$f + 9" | bc) ; 
sed -n "8,$n p" CONTCAR  >> Cluster/POSCAR 

cd Cluster
pwd
sbatch run.sh
ls
cd ..
