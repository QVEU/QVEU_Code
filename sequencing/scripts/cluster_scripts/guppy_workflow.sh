# Benjamin Adam Catching
# 2022-12-14
# NIH-NIAID-LVD
# Quantitative Virology and Evolution Unit


# Make sure you are in the correct directory relative to the files
# Works so far with all these lines typed in command line, untested with submitting bash script

# Log onto development node with 1 gpu
qrsh -l gpu=1,v100
# Load guppy with gpu
module load guppy_gpu/6.0.6
# Run guppy basecaller from command line, make sure the flowcell and kit combination are useable by typing -> guppy_basecaller --print_workflows
guppy_basecaller -i fast5_pass/barcode13 -s guppy_basecalled/barcode13 --flowcell FLO-MIN106 --kit SQK-NBD110-24 -x cuda:0
# -i is the folder where the fast5 files are
# -s is the folder where the fastq files should be saved (make this file before running)
# -x for using the gpu (cuda:0 for 1 gpu, cuda:0:1 for 2 gpus etc.)
