#!/bin/bash -e

#SBATCH --job-name      RepRunning
#SBATCH --time          3:00:00
#SBATCH --ntasks        1
#SBATCH --cpus-per-task 30
#SBATCH --array         1-20
#SBATCH --mem-per-cpu   4GB
#SBATCH --partition=####
#SBATCH --output        RepRunning.%j.out # Include the job ID in the names of
#SBATCH --error         RepRunning.%j.err # the output and error files
#SBATCH -A #### 
cd /####/sims_complete_rerun_31oct
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
Rscript replicate_running_nonmpi_multigen_tm1.R
#Rscript replicate_running_nonmpi_multigen.R
#Rscript replicate_data_sorting.R
echo "R finished."


# Define the base paths and file patterns
input_data_base="tm_1_2Nov_results_multigen_replicate_"
envdata_base="rerun_1nov_rerun_sec_varyring_envs_multigen_replicate_"
output_base="tm_1_rerun_2nov_largererun_data_analysed_replicate_"

# Get the replicate number from the SLURM_ARRAY_TASK_ID
replicate_num=${SLURM_ARRAY_TASK_ID}

# Construct the file paths
input_file="${input_data_base}${replicate_num}.RDS"
envdata_file="${envdata_base}${replicate_num}.RDS"
output_file="${output_base}${replicate_num}.RDS"
echo $input_file
echo $envdata_file
echo $output_file

# Run the R script with the specified files
Rscript ./generate_analysed_RDS_updated.R $input_file $envdata_file $output_file

