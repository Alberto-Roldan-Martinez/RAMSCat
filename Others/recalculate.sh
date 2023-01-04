#!/bin/bash

rm Measured_Relaxed.dat  Measured_Static.dat  Predicted.dat
head -30 Au2/1/Predicted.dat >> Predicted.dat;
cp Predicted.dat Measured_Static.dat
#cp Predicted.dat Measured_Relaxed.dat

for i in 2 3 4 5 6 7 8 9 10 20
 do cd Au$i;
  re='^[0-9]+$'; 
  for j in $(ls); 
    do if [[ $j =~ $re ]] 
      then cd $j;
      pwd
      rm CONTCAR.vasp  Optimisation.txt  output  POSCAR  Predicted.dat  RAMSCat.out 
      cp Validation/Static/POSCAR .
      python3 ~/Software/RAMSCat/Predicting.py Au POSCAR MgO 8 8 2 $i >> output
      more output
      tail -1 Predicted.dat >> ../../Predicted.dat;
      tail -1 Validation/Static/Measured.dat >> ../../Measured_Static.dat ;
#      tail -1 Validation/Relaxed/Measured.dat >> ../../Measured_Relaxed.dat ;
      cd ..;
     fi ;
    done;
  echo "" >> ../Predicted.dat;
  echo "" >> ../Measured_Static.dat ;
#  echo "" >> ../Measured_Relaxed.dat;
 cd ..;
done;
python3 ~/Software/RAMSCat/Pre_Data_Collection/Plot.py Measured_Static.dat Predicted.dat
