#!/bin/bash


for i in Au2 Au3  Au4  Au5  Au6  Au7  Au8  Au9 Au10
 do cd $i; 
  re='^[0-9]+$'; 
  for j in $(ls); 
    do if [[ $j =~ $re ]] 
      then cd $j;
      pwd
      rm output Predicted.dat RAMSCat.out
      python3 ~/Software/RAMSCat/Predicting.py Au CONTCAR.vasp MgO 8 8 2 2 >> output
      cd ..;
     fi ;
    done;
 cd ..;
done;
~/Software/RAMSCat/Others/cross_relation.sh
python3 ~/Software/RAMSCat/Pre_Data_Collection/Plot.py Measured_Relaxed.dat Predicted.dat