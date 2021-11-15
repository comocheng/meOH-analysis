# script to do the 5000 RMG runs on discovery
# Sevy Harris 9/8/2021

import numpy as np
import job_manager
import os
import glob


# WARNING - this will fail if M%N != 0


skip_completed_runs = False  # set to false to overwrite RMG runs that completed

working_dir = "/scratch/westgroup/methanol/perturb_5000/"
# working_dir = "/home/moon/rmg/fake_rmg_runs/"
if not os.path.exists(working_dir):
    os.mkdir(working_dir)

M = 5000  # total number of times to run RMG
N = 50  # number of jobs to run at a time

for i in range(0, M, N):
    sbatch_index = int(i / N)
    range_max = np.amin([i + N, M])
    last_index = range_max - 1
    job_indices = [a for a in range(i, range_max)]

    # max slurm array index is 1000, so after that, subtract multiples of 1000
    task_id_offset = int(i/1000) * 1000
    
    # Write the job file
    fname = f'ct_graaf_runs_{i}-{last_index}.sh'
    script_folder = "ct_run_scripts"
    jobfile = job_manager.SlurmJobFile(full_path=os.path.join(working_dir, script_folder, fname))
    jobfile.settings['--array'] = f'{i - task_id_offset}-{last_index - task_id_offset}'
    jobfile.settings['--job-name'] = fname
    jobfile.settings['--error'] = os.path.join(working_dir, script_folder, f'ct_graaf_job_error{sbatch_index}.log')
    jobfile.settings['--output'] = os.path.join(working_dir, script_folder, f'ct_graaf_job_output{sbatch_index}.log')
    jobfile.settings['--mem'] = f'20Gb'
    jobfile.settings['--cpus-per-task'] = '4'
    
    content = ['# Define useful bash variables\n']

    
    content.append(f'SLURM_TASK_ID_OFFSET={task_id_offset}\n')
    content.append('RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID + $SLURM_TASK_ID_OFFSET)))\n')
    rmg_run_dir = os.path.join(working_dir, "run_${RUN_i}")
    content.append(f'CSV_FILE="{rmg_run_dir}/cantera/ct_graaf_analysis.csv"\n')
    content.append(f'CT_FILE="{rmg_run_dir}/cantera/chem_annotated.cti"\n')
    
    # skip if csv file already exists
    if skip_completed_runs:
        content.append('if test -f "$CSV_FILE"; then\n')
        content.append('echo "skipping completed run ${RUN_i}"; exit 0\n')
        content.append('fi\n\n')
        
    
    # run the analysis script
    content.append('# remove the old CSV file\n') 
    content.append('rm -f $CSV_FILE\n') 
    content.append('# Run the Cantera analysis\n') 
    content.append(f'python "/scratch/westgroup/methanol/meOH-analysis/cantera_simulations/uncertainty/analyze_graaf_runs.py" $CT_FILE\n')
    jobfile.content = content
    jobfile.write_file()
   
 
