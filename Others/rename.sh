#!/bin/bash

tar -xvf Structures.tar
 rm Structures.tar 
 cd Structures 

re='^[0-9]+$'
 for i in $(ls); do 
	if [[ $i =~ $re ]] ; then 
		a=$(($i + 1))
		 mv $i a$a
	fi
 done;
 for i in $(ls); do
	mv $i ${i:1}
 done 
ls
