#!/bin/bash


rm Measured_Relaxed.dat  Measured_Static.dat  Predicted.dat
head -30 Au2/1/Predicted.dat >> Predicted.dat;
cp Predicted.dat Measured_Static.dat
cp Predicted.dat Measured_Relaxed.dat

for i in Au2 Au3  Au4  Au5  Au6  Au7  Au8  Au9 Au10
 do cd $i; 
  re='^[0-9]+$'; 
  for j in $(ls); 
    do if [[ $j =~ $re ]] 
      then cd $j;
#      pwd
      tail -1 Predicted.dat >> ../../Predicted.dat;
      tail -1 Validation/Static/Measured.dat >> ../../Measured_Static.dat ;
      tail -1 Validation/Relaxed/Measured.dat >> ../../Measured_Relaxed.dat ;
      cd ..; 
     fi ;
    done;
  echo "" >> ../Predicted.dat;
  echo "" >> ../Measured_Static.dat ;
  echo "" >> ../Measured_Relaxed.dat;
 cd ..;
done;
