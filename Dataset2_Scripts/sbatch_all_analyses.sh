#!/bin/bash
#Submit jobs for all differential abundance methods

for dir in *
do
	if [ -d "$dir" ]; then
		cd $dir
		for job_file in *.job
		do
			echo "Submitting job file $job_file"
			sbatch $job_file
			sleep 15s
		done
		cd ..
	fi
done
