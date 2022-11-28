## Cell Ranger Single Cell Pipeline

1. Make fastQ files using the BCLs from Illumina sequencing;
2. Make Cell Ranger counts using the Fastqs from step 1;
3. Make Cell Ranger aggregate using the counts output of the appropriate conditions. 

### 1 - Make fastQs with mkfastq Cell Ranger

Use shell script as example to set up yours on QVEU_Code/SingleCell/CellRangerScripts folder + the pathway to the BCL file (the pathway will goes just until the name of the flowcell - the software knows how to get from there). 

example: 

    qsub mkfastq_QVEU0034.sh /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0034_NextSeq2K_SpikeInTest2_withadaptors_3/221123_VH01023_12_AAANG3NHV
    
Note: you can also use **sh** but then you can't disconect your computer from the internet (or even which sources - cable vs wi-fi) otherwise you will terminate your job.

