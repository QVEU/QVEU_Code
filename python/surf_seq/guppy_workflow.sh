# Benjamin Adam Catching
# 2023-03-01
# NIH-NIAID-LVD
# Quantitative Virology and Evolution Unit

"""
This pipeline is for automatically taking a Nanopore run with the SQK-NBD112-96 kit and 
returning 
"""
# Make sure you are in the correct directory relative to the files
# Works so far with all these lines typed in command line, untested with submitting bash script

# Log onto development node with 1 gpu
qrsh -l gpu=1,v100
# Load guppy with gpu
module load guppy_gpu/6.0.6

# Make appropriate folders
touch barcoding
touch guppy_basecalled 
touch Q20
touch sam
touch duplex_calls

# Run guppy basecaller from command line, make sure the flowcell and kit combination are 
# useable by typing -> guppy_basecaller --print_workflows
guppy_basecaller -i fast5 -s guppy_basecalled \
--flowcell FLO-MIN106 --kit SQK-NBD112-96 -x cuda:0
# -i is the folder where the fast5 files are
# -s is the folder where the fastq files should be saved (make this file before running)
# -x for using the gpu (cuda:0 for 1 gpu, cuda:0:1 for 2 gpus etc.)

# Barcoding can be done from the output folder, folders for each barcode 
# vwill be generated within the barcoding folder
guppy_barcoder -i guppy_basecalled/pass -s barcoding -x cuda:0

# Consolidate all fastq files into one fastq file for a single mapping to .sam files
python fastq_consolidate.py 

# Map reads to .sam files
bash minimap_96_sample_loop.sh

# Map reads to mutation observances
bash loop_96_sam_to_Q20.sh

# Experimental section, please don't run quite yet
#duplex_tools pairs_from_summary <sequencing_summary.txt/dorado.bam> <output directory>

#duplex_tools filter_pairs <pair_ids.txt> <fastq directory/dorado.bam>

#guppy_basecaller_duplex -i . -r -s duplex_calls -x 'cuda:0' --chunks_per_runner 416 \
# --duplex_pairing_mode from_pair_list --duplex_pairing_file pair_ids_filtered.txt
