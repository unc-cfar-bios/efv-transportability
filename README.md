# efv-transportability

#Transportability analysis with multiple imputation and percentile-based bootstrap 95% CI:

Each run of the R script bootMI_4cluster creates B=5 bootstrap iterations. 
The shell script contains a command to run a Slurm job array of 20 jobs (SBATCH --array=1-20), 
with each run of the shell script producing 100 (=5*20) bootstraps. The shell script was run 
repeatedly in the console (100 times with 2 second delays in between) to produce 10,000 bootstrap resamples 
in an appended CSV file.


Instructions for Installing R Packages on a Linux Cluster
STEP 1: Load the R module 
Note: R may be capitalized or lowercase depending on your cluster’s naming convention.
module load R    or    module load r

STEP 2: Start a new R session
R 
STEP 3: Install “batch” package
install.packages(“batch”) 
Note: If there is a warning similar to 
Warning in install.packages("batch") :
'lib = "/nas/longleaf/apps/r/3.3.1/lib64/R/library"' is not writable

then follow the prompts to use a personal library instead, and select the applicable CRAN mirror (select the mirror nearest to your current location).
STEP 4: Load and attached the “batch” package
library(batch)

STEP 5: Install the necessary package(s)
Example: install.packages(“boot”) 
STEP 6: End the R session
quit()





