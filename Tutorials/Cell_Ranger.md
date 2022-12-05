## Cell Ranger Single Cell Pipeline

1. Make fastQ files using the BCLs from Illumina sequencing;
2. Make Cell Ranger counts using the Fastqs from step 1;
3. Make Cell Ranger aggregate using the counts output of the appropriate conditions. 

### 1 - Make fastQs with mkfastq Cell Ranger

Use shell script as example to set up yours on QVEU_Code/SingleCell/CellRangerScripts folder + the pathway to the BCL file (the pathway will goes just until the name of the flowcell - the software knows how to get from there). 

example: 

    qsub mkfastq_QVEU0034.sh /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0034_NextSeq2K_SpikeInTest2_withadaptors_3/221123_VH01023_12_AAANG3NHV
    
Note: you can also use **sh** but then you can't disconect your computer from the internet (or even which sources - cable vs wi-fi) otherwise you will terminate your job.

### 2 - Make Cell Ranger counts

Use shell script as example to set up yours on QVEU_Code/SingleCell/CellRangerScripts folder. Ideally, you will add the pathway to the folder you want to save your counts or run it from the folder that it should be saved.
*fqpath - you will add the pathway for all fastqs you have for this sample. You can add multiple fastqs from multiple runs. You should add the fastq pathway until the folder that has the number of the flow cell, if you have more the more than one run add both pathways separated by comma (,).
*template - add the pathway to the Reference genome to your sample (remember if you have a custom reference here!) 

example of the shell script: 

    #!/bin/sh
    #$ -pe threaded 16
    #$ -l h_vmem=8G
    #$ -m bea
    #$ -M silvaaln@nih.gov

    #USAGE: for i in *fastq; do; qsub CountJob.sh $i AAANG5WHV_12; done;

    module load cellranger/7.0.0
    cd /hpcdata/lvd_qve/Projects/QVEU_Seq_0025_0031_0034_Analysis/
    fqpath=/hpcdata/lvd_qve/Projects/QVEU_Seq_0025_Analysis/AAANG5WHV_12/,/hpcdata/lvd_qve/Projects/QVEU_Seq_0025and0031_Analysis/AAANGTJHV_12/,/hpcdata/lvd_qve/Projects/QVEU_Seq_0025_0031_0034_Analysis/AAANG3NHV_12
    template=/hpcdata/lvd_qve/QVEU_Code/sequencing/template_fastas/REF_human-NPEVs/hg38_CVB3-Nancy/

    runID=${1}
    echo running CVB3.
    cellranger count --sample=CVB3_1 --transcriptome=$template --id=CVB3_1 --fastqs=${fqpath}${runID}/outs/fastq_path/ --localcores=16 --localmem=125

You can use the sh (interactive) or qsub (unsupervised job) to run the shell script - it takes more than 2 hour to run (depending on the size of your sequecing runs - in here I used 3 sequecing runs, from NextSeq P3 flow cells).

### 3 - Make Cell Ranger aggregate using the counts output of the appropriate conditions

Use shell script as example to set up yours on QVEU_Code/SingleCell/CellRangerScripts folder.

{1} is your .csv file

example of csv file: 

![image](https://user-images.githubusercontent.com/97693929/205693790-6f5aeeab-2ddb-40f4-8265-361541d6a541.png)


{2} is the name you want to give to that folder

example:

    sh AggrJob.sh aggr_CVB3.csv CVB3

You can use qsub as well for a non supervised job. 

### Exploratory analysis in the web_summary and Loupe browser

**1. Web_summary results** of each sample will give a couple of important parameters for you go foward (or not) with you analysis

*Summary tab:

  1.1 - **Estimated Number of Cells** (this number will be around of the target cells you choose during cells preparation)
  1.2 - **Mean Reads per Cell** (which should be above 20k per cell - at least), you should check the requeriments according to which chemestry you are using (3' or 5'). Sometimes sequecing more will provide better covarege per cell = more confidence of the transcriptome you are detecting. For 5' kit, if you are not using VDJ sequecing, just gene expression, there is a note for 50k per cell on the end of the 10X protocol. 
  1.3 - **Median Genes per Cell** (you can detect between 2k to 7.5k genes in one cell)
  1.4 - **Sequecing Saturation** - this value should be around 100%. If you are not it could be 2 things - you need to sequence more to be able to saturate the heterogeneity in your sample. However, if you are working with low complexity cells (such as monocultered cells) this could be challenging. One good parameter that can help to decide if you need more sequecing is in the gene expression tab (Median Genes per Cell), which will show if you are reaching a plato of genes detected per cell. 
    
*Gene expression tab:
 
  1.5 - **tSNE of UMIs counts** UMI means number of strands those transcript where captured, if you have low number of genes per cell or/and low number UMI per cell it could mean that "cell" is not trully a cell but just background. 
    ![newplot](https://user-images.githubusercontent.com/97693929/205698722-ec53945d-dd56-40ff-9c13-253dec20fa07.png)

  1.6 - **tSNE of cells by clustering** This could be usefull for further analysis giving you a clue of the number of expected cluster in your sample (however this clusters will need validation according to the genes markers the exhibit - be aware!)
  ![newplot (1)](https://user-images.githubusercontent.com/97693929/205700478-d76fe044-79c4-44fa-b70a-7e0edba2f5b3.png)
  
  1.7 - **Top Features by Cluster (Log2 fold-change, p-value)** - If you have prior knowledge of your sample, you may search if the expected genes pop up on your clusters, or if you don't, you can use this list to check about their "functions" briefly. 

**2: Loupe browser results**

You will have cloupe.cloupe results from step 2 and step 3 of this pipeline, from each individual sample or from the aggregate of each appropriate condition (Mock and infection, for example) to explore.

Using the aggreg results, for example, you may run the globally distinguishing calculation between the different sample IDs. This will generate a list of up-regulated genes per sample, which works as a bulk RNAseq. Those genes would be good to validate your sequencing results, since you will probably find a bulk RNAseq in the literature. 

Check for "viral reads" on your MOCK sample, if you find some counts - it could indicate a "background" (human reads that are align wrong to the viral gene - however, if the counts are high be careful with "contamination") to you "regress out" of your infected condition. Otherwise, you could be superestimating the true viral reads of your infected condition. 


Note: all this data is for a brief check exploration, you will use them to decide if you already achieve the sequecing depth you want or you need to sequence more your samples. Remember all this data was performed by unsupervised setting of Cell Ranger, they are insuficient for publications, for example. You will need another single-cell pipeline for future publications such as Seurat, for example.
