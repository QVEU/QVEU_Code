BioWULF Tutorial or: How I Learned to Stop Worrying and Love the Code
-----------------------------------------------------------------------
_Jack Dorman_
7/21/22

BioWULF is the NIH's main High Performance Computing cluster. The following tutorial is based on the BioWULF User Guide found at https://hpc.nih.gov/docs/userguide.html as well as example scripts provided by Nidia plus a little bit of trial and error. The BioWULF website also has some tutorial videos if that's more your style. You should come into this with some familiarity with Linux commands. In addition, it uses slurm for job submission which means it has identical commands to other computer clusters that use slurm.

## Getting in and Setting Up
### Log onto BioWULF
Use ssh command to remotely log onto the BioWULF HPC. The first time you log in,
but not again in subsequent log ins, you'll be asked to confirm. Enter yes then you should be prompted for your NIH passphrase.
```{bash}
ssh username@biowulf.nih.gov
```

### Moving Files to and From BioWULF
Use scp to transfer files from your local computer to BioWULF.
```{bash}
scp ~/local/path/Job_Script.sh biowulf.nih.gov:/home/username/ #from local path to BioWULF
scp biowulf.nih.gov:/home/username/results.fasta ~/local/path/ #from BioWULF to local
```
To transfer a whole directory with all its contents use -r.
```{bash}
scp -r ~/local/path/Job_Script.sh biowulf.nih.gov:/home/username/
```


## Running Jobs
### Script Set Up
.sh scripts for running jobs should have the following at the top. Other parameters can be included but these should be the minimum. See example scripts for commands for different jobs.
```{bash}
#!/bin/bash
# this file is Job_Script.sh
#SBATCH --job-name=*Short Name*
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="email.address@nih.gov" #emails you when the job is submitted, when it ends, or if it fails

set -e #ends job if a command or signal fails

```

### Submitting jobs
Use sbatch to submit a job. The minimum command is as follows. Upon submission, a file called slurm-jobID.out will be created in the directory where sbatch is run. It will contain any standard output or standard error messages.
```{bash}
sbatch ~/path/to/script/job_Script.sh #if the script is in the same directory that you're running the sbatch then just enter the script name
```
Other parameters can be entered and the full list is in the BioWULF user guide. The following sets the run time at a maximum of 3 days, specifies the type of node to be used, and specifies how many cpus to use for each task. When the job is first submitted it is given a job ID number and it will likely wait in the queue before running. 
```{bash}
sbatch --time=3-00:00:00 --constraint=ibhdr100 --cpus-per-task=10 ~/path/to/script/Job_Script.sh
```
To check on the job use the following. You will see, among other things, the job ID number, the name, if it's running (R) or pending (PD), and how long it's been running.
```{bash}
squeue -u *username*
```
To cancel a job:
```{bash}
scancel *jobID*
```


## R
### Call an interactive R session
For running R commands interactively such as to install packages in personal library.
```{bash}
sinteractive #this will take a moment to load
module load R/4.2 #this is the latest installed version of R but other versions can be called
R
```

### Install Packages in Personal Library
Example below is installing DECIPHER used in codon alignment. Afterwards you will have a new directory in your home directory unless another location is specificed. The directory will be ~/R/_version_/library
```{r}
BiocManager::install("DECIPHER")
```

### Running a Job with Rscript
In addition to the Rscript, which is made like any other Rscript, a separate .sh script must be made to run it. Below is an example script which runs the Rscript called CodonAlign_Decipher.R. This .sh script should be submitted like any other.
```{bash}
#!/bin/bash
# this file is rJob.sh
#SBATCH --job-name=*Short Name*
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="email.address@nih.gov" 

set -e 

module load R/4.2.0
Rscript ~/CodonAlign_Decipher.R
```

## Troubleshooting
If the run fails, first check the slurm.out file. If no clear error is shown, go to the User Dashboard on the BioWULF website then Job Info and finally select the run ID. The most likely cause is running out of memory. On the cluster, the following command can also give more thorough info on the run:
```{bash}
jobhist *JobID*
```
