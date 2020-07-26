#!/bin/sh

#SBATCH --job-name=coolest_job_ever
#SBATCH --time=11-00:00:00
#SBATCH --mail-user=your_email@domain.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p general
#SBATCH --cpus-per-task=1
#SBATCH --array=1-20

R --vanilla </PATH_HERE/bootMI_4cluster.R
