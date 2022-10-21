#!/bin/bash

# validation n_atoms


re='^[0-9]+$'; 
for j in $(ls); 
	do if [[ $j =~ $re ]] 
		then cd $j;
		if [ ! -d "Validation" ];
			then cp -rf ../../Au2/1/Validation/ .
			cp CONTCAR.vasp Validation/Relaxed/POSCAR
                        cp CONTCAR.vasp Validation/Static/POSCAR
			cd Validation/Relaxed/
			rm CONTCAR O* Measured.dat
			pwd
			sbatch run.sh

			cd ../Static/
                        rm O* CONTCAR Measured.dat
			pwd
			sbatch run.sh

			if [ -f "Cluster/POSCAR" ];
			        then rm Cluster/POSCAR  Cluster/O* Cluster/CONTCAR
			fi
			sed -n "1,5 p" POSCAR >> Cluster/POSCAR
			sed -n "6 p" POSCAR | awk '{print $1}' >> Cluster/POSCAR
			f=$(sed -n "7 p" POSCAR | awk '{print $1}') # >> Cluster/POSCAR   
			echo $f >> Cluster/POSCAR
			n=$(echo "$f + 9" | bc) ;
			sed -n "8,$n p" POSCAR  >> Cluster/POSCAR

			cd Cluster
			pwd
			sbatch run.sh
			cd ../../../
		fi
		cd ../
	fi
done
ls

