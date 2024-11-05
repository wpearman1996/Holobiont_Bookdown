#!/bin/bash -e

#SBATCH --job-name      RepRunning
#SBATCH --time          1:00:00
#SBATCH --ntasks        1
#SBATCH --cpus-per-task 1
#SBATCH --array         1-20
#SBATCH --mem-per-cpu   8GB
#SBATCH --output        RepRunning.%j.out # Include the job ID in the names of
#SBATCH --error         RepRunning.%j.err # the output and error files
#SBATCH -A #### 

cd /#####/holobiont_first_paper_sim_reps/sims_28Aug24/evi_analysis

# Load the R module
module load R

# Define the base paths and file patterns
input_data_base="evi_1gen_wn_largererun_results_multigen_replicate_4Nov_"
envdata_base="rerun_1nov_rerun_sec_varyring_envs_multigen_replicate_"
output_base="evi_1gen_wn_largererun_data_analysed_replicate_4Nov"

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
Rscript generate_analysed_RDS.R $input_file $envdata_file $output_file
