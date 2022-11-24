#!/bin/bash


re='^[0-9]+$'
declare -i n=1
#n=1
for i in $(ls); do
	if [[ $i =~ $re ]]; then
		n=1
		while [[ -d $n ]]; do
#			n=$((n+1))
			n+=1
		done
		mv $i $n
	fi
done;
ls
