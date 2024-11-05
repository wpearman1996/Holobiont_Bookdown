#!/bin/bash -e

#SBATCH --job-name      RepRunning
#SBATCH --time          02:00:00
#SBATCH --ntasks        1
#SBATCH --cpus-per-task 15
#SBATCH --array         1-20
#SBATCH --mem-per-cpu   2GB
#SBATCH --partition=####
#SBATCH --output        RepRunning.%j.out # Include the job ID in the names of
#SBATCH --error         RepRunning.%j.err # the output and error files
#SBATCH -A #### 
cd /#####/holobiont_first_paper_sim_reps/sims_28Aug24/evi_analysis
module load R/4.2.1-gimkl-2022a
# need MPI
#module load gimkl/2022a
module load GCC

# Help R to flush errors and show overall job progress by printing
# "executing" and "finished" statements.
echo "Executing R ..."
# Our R has a patched copy of the snow library so that there is no need to use
# RMPISNOW.

#srun Rscript replicate_running.R
#Rscript replicate_running_nonmpi.R
Rscript replicate_running_nonmpi_multigen.R
#Rscript replicate_data_sorting.R
echo "R finished."
