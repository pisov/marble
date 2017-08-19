#!/bin/bash

#SBATCH -p long.p               # partition (queue)
#SBATCH --job-name="produc"          # name
#SBATCH -N 8                      # number of nodes
#SBATCH -n 256                    # number of cores
#SBATCH -t 5-0:00                 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out        # STDOUT
#SBATCH -e slurm.%N.%j.err        # STDERR

module load openmpi

folder=0
for a in {500..1000..25}
do	
	folder=$(($folder+1))
	var=$(awk 'BEGIN { printf "%0.5f", "'$a'"/100000}')
#	rm -r $a

	mkdir -p $folder
#	echo $folder
	cp `ls -p | grep -v / `  $folder/
	cd $folder
#        sed   "/^pbvac =/s/=.*/=$var/" marble.in2 
        sed -i  "/^pbvac =/s/=.*/=$var/" marble.in 
        mpirun -np 256 marble.x >> run.dat
	rm $(ls -I "*.dat" -I "*.in" -I "*chkp*")

	cd ..
done
