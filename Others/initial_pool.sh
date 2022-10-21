#!/bin/bash


element=Au

for i in 10 # 2 3 4 5 6 7 8 9 10; 
	do #cd Au$i;
	a=$(grep energy 1/RAMSCat.out | awk '{print $3}'); 
	s=$(grep energy 1/RAMSCat.out | awk '{print $6}');
	echo "$i" >> pool.dat;
	echo "Energy = $a Dir = 1 Sphericity = $s" >> pool.dat
	j=1
	while [ $j -le $i ]; do
		f=$(echo "$j + 3" | bc) ;
		position=$(sed -n "$f p" 1/RAMSCat.out) 
		echo "$element $position" >> pool.dat
		((j++))
	done
	cd ..
done

