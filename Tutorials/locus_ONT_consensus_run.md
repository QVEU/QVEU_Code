# Generating consensus sequences from ONT runs. 
## Author: Patrick T. Dolan, Ph.D.
## Date: 5/6/22

The read directories on the nanopore computer look like: 

<img width="384" alt="Screen Shot 2022-05-06 at 8 02 54 AM" src="https://user-images.githubusercontent.com/10180619/167127791-f40fb881-e4d2-400a-988b-8ee3a133f1df.png">

Move these files to the _locus.niaid.nih_ cluster [see tutorial here.](https://github.com/QVEU/QVEU_Code/blob/main/Tutorials/locus_tutorial.md#3-transferring-files-to-the-cluster-with-rsync)

### 1. In your home directory clone the `QVEU_Code` GitHub repo. 
```
git clone https://github.com/QVEU/QVEU_Code
```

### 2. Edit the run file, `~/QVEU_Code/sequencing/scripts/cluster_scripts/nanopore_consensus_run_<yourinfo>.sh`:

The run file is just a simple bash script that contains all the necessary modules and the script. 

``` bash
#!/bin/bash

#nanopore_consensus_run_<yourinfo>.sh
#edit this script for individual runs on the cluster.
#open an interactive job with 'qrsh' and then sun `sh minion_run_<yourinfo>.sh`. 

#load required modules...
module load minimap2
module load samtools

#run minimap and samtools to generate consensus
sh ~/QVEU_Code/sequencing/scripts/nanopore_consensus_scripts/minimap_minion_consensus.sh ~/QVEU_Code/sequencing/template_fastas/<yourtemplate>.fasta ~/<path to your data>/fastq_pass/
```

Running this file `sh ~/QVEU_Code/sequencing/scripts/cluster_scripts/nanopore_consensus_run_<yourinfo>.sh` will set up the environment on locus with the necessary modules, and then execute the minimap mapping script, `minimap_minion_consensus.sh`, on the defined template and fastq directory. 

### 3. Navigate from the top directory to the 'fastq_pass' directory on the computer you will run your mapping script:
``` bash
cd <my_lib_name>/<my_sample_name>/<run_id>/fastq_pass/
```



``` bash
cd <my_lib_name>/<my_sample_name>/<run_id>/fastq_pass/
```

