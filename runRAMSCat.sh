#!/bin/bash
#SBATCH --account=scw1039
#SBATCH -p c_compute_chemy1		# -----if it is more than 40 tast, it should be in compute
#SBATCH -o file.%J
#SBATCH -e error.%J
#SBATCH --ntasks-per-node=40      # tasks to run per node
#SBATCH -n 40                     # number of parallel processes (tasks)
#SBATCH -t 72:00:00                # time limit
#SBATCH -J np2                # Job name


# Load the environment
   	source ~/venv/bin/activate
	module load python/3.7.0
	module load compiler/gnu/7/3.0
	module load ase/3.20.1

# Setting directories and defining variables
        export PYTHONPATH=$PYTHONPATH:"/home/c.sacar6/bin/RAMSCat/"
        NNODES=$SLURM_NNODES
        NCPUS=$SLURM_NTASKS
        PPN=$SLURM_NTASKS_PER_NODE
        MYPATH=$SLURM_SUBMIT_DIR
        WDPATH=/scratch/$USER/$SLURM_JOBID
# Setting environment and creating scratch directory
        mkdir -p ${WDPATH} ; cd ${WDPATH} ;
        cp -rf ${MYPATH}/* .
        echo ${MYPATH} >> output

# Launch the parallel job Using ncpus processes
        env; echo RAMSCat Start Time is `date` running NCPUs=$NCPUS PPN=$PPN
        start="$(date +%s)"

    python3 runPredicting.py >> output 

        echo RAMSCat Finish Time is `date` ; stop="$(date +%s)" ; finish=$(( $stop-$start ))
        echo RAMSCat $SLURM_JOBID Job-Time  $finish seconds
	deactivate ## deactivates the environment

# Copy output data to home
	mkdir Structures; mv 1* 2* 3* 4* 5* 6* 7* 8* 9* Structures;
	for i in $(grep Energy RAMSCat_Summary.txt | awk '{print $6}'); do cp -rf Structures/$i .; done;
	tar -zcf Structures.tar Structures; rm -rf Structures
        mv ${WDPATH}/* ${MYPATH}/
	cd ${MYPATH}/
	echo "seff $SLURM_JOBID | grep \"Efficiency\" >> RAMSCat_Summary.txt" >> efficiency.sh ; chmod +x efficiency.sh
        rm -rf ${WDPATH}

# To remove all the NUMERIC folders
#	re='^[0-9]+$'; for i in $(ls); do if [[ $i =~ $re ]] ; then rm -rf $i; fi; done
