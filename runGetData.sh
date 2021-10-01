#/bin/sh


  dir=$(pwd)

#	file=$(find . -name "x*" -type d)
	for i in $(echo "x*");
       	do 
		cd $i
		echo $(pwd)
		python3 /home/alberto/Software/OTHER/NeuralNetwork/GetData.py $dir
		if [ -f "data.dat" ]; then
			cat data.dat >> ../Data.tmp
			rm data.dat
			mv labels.txt ../.
		fi
		cd ..
       	done
	cat labels.txt Data.tmp >> Data.dat
	rm labels.txt Data.tmp
#	python3 /home/alberto/Software/OTHER/NeuralNetwork/Plot.py Data.dat



