#/bin/sh


  dir=$(pwd)

	file=$(find . -name "x*" -type d)
	for i in $file;
       	do 
		cd $i
		python3 /home/alberto/Software/OTHER/NeuralNetwork/GetData.py $dir
		if [ -f "data.dat" ]; then
			cat data.dat >> ../Data.tmp
			rm data.dat
			mv labels.txt ../.
      python3 ~/Software/OTHER/NeuralNetwork/CatStructure.py
		fi
		cd ..
       	done
	cat labels.txt Data.tmp >> Data.dat
	rm labels.txt Data.tmp
	python3 /home/alberto/Software/OTHER/NeuralNetwork/Plot.py Data.dat

# to normalise -- implemented before Ebinding
#	python3 /home/alberto/software/OTHER/NeuralNetwork/ENormalisation.py Data.dat
#	python3 /home/alberto/software/OTHER/NeuralNetwork/Plot.py NData.dat


