# efv-transportability

#Transportability analysis with multiple imputation and percentile-based bootstrap 95% CI:

Each run of the R script bootMI_4cluster creates B=5 bootstrap iterations. 
The shell script contains a command to run a Slurm job array of 20 jobs (SBATCH --array=1-20), 
with each run of the shell script producing 100 (=5*20) bootstraps. The shell script was run 
repeatedly in the console (100 times with 2 second delays in between) to produce 10,000 bootstrap resamples 
in an appended CSV file.



