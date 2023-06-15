# Benjamin Adam Catching
# 2023-03-02
# NIH-NIAID-LVD
# Quantitative Virology and Evolution Unit

"""
Leptocottus armatus (Pacific staghorn sculpin), v0.2
From aligned, phred-score = 20 cutoff, Q20.txt files, call each consensus
mutation in each replicate of a 8-replicate, 12-condition passaging experiment
"""

# Import modules
import numpy as np
import pandas as pd
import glob
from Bio import Seq

print('Reading in files')
new_P1_dir = 'data/sequences/2023-02-06_passage_1/Q20/*.txt'
new_P2_dir = 'data/sequences/2023-02-27_passage_2/Q20/*.txt'
new_P3_dir = 'data/sequences/2023-03-14_passage_3/Q20/*.txt'
new_P4_dir = 'data/sequences/2023-04-06_passage_4/Q20/*.txt'
new_P5_dir = 'data/sequences/2023-05-30_passage_5/Q20/*.txt'

directories = [new_P1_dir, new_P2_dir, new_P3_dir, new_P4_dir,
               new_P5_dir]

# Iterate through Q20 directories
passage_files = [sorted(glob.glob(directory)) for directory in directories]

# List comprehension of the data
Q20_df =[[pd.read_csv(x, delimiter='\t') for x in y] for y in passage_files]

# Define conditions
strains = ['US/MO/14-18947', 'Fermon']
cell_types = ['SH-SY5Y', 'RD', 'A549']
temps = ['33', '37']
print('Converting to DataFrame')
for j, passage in enumerate(Q20_df):
    for i, replicate in enumerate(passage):
        batch = i // 8
        replicate_num = i % 8
        replicate_chr = chr(ord('@') + replicate_num+1)
        strain_num = batch % 2
        strain = strains[strain_num]
        cell_num = batch // 4
        cell_type = cell_types[cell_num]
        temperature_num = batch // 2 % 2
        temperature = temps[temperature_num]
        
        replicate['passage'] = str(j+1)
        replicate['cell'] = cell_type
        replicate['strain'] = strain
        replicate['temperature'] = temperature
        replicate['replicate'] = replicate_chr
        replicate['nucleotide position'] = [x + 1 for x in list(replicate.index)]

# Concatanate Q20 file
Q20_df_list = pd.concat([pd.concat(Q20_df[i]) for i in range(len(Q20_df))])
Q20_df_list = Q20_df_list.astype({"passage": int, "cell": str, "strain": str, 
	"replicate": str, "temperature":int})
Q20_df_list = Q20_df_list.rename(columns={'0': 'A', '0.1': 'C', '0.2': 'G', 
	'0.3': 'T', '0.4': 'N'})

Q20_df_list = Q20_df_list.drop(['1','12','14','8','2','4','11','3','13'], axis=1)
Q20_df_list['coverage'] = Q20_df_list['A'] + Q20_df_list['C'] +  Q20_df_list['G'] + Q20_df_list['T']

#Q20_df_list.to_csv('data/EV-D68_all_passage_data_v3.csv')
print('Melting DataFrame into nucleotide-position')
melted_df = pd.concat([Q20_df_list[[x, 'passage', 'cell', 'strain', 'temperature',
         'replicate', 'nucleotide position', 'coverage']].melt(['passage', 'cell', 'strain', 'temperature',
         'replicate', 'nucleotide position', 'coverage']) for x in ['A', 'C', 'G', 'T']])

melted_df['percent'] = melted_df['value'] / melted_df['coverage']

# Need a list of Fermon protein start and stops
fermon_protein_pos = {"5' UTR":0, 'VP4':731, 'VP2':938, 'VP3':1682, 'VP1':2387,
                     '2A':3314, '2B':3755, '2C':4052,
                     '3A':5042, '3B':5309, '3C':5375, '3D':5924,
                     "3' UTR":7298}

# Need a list of US/MO/14-18947 protein start and stops
MO_protein_pos = {"5' UTR":0, 'VP4':697, 'VP2':904, 'VP3':1648, 'VP1':2353,
                     '2A':3280, '2B':3721, '2C':4018,
                     '3A':5008, '3B':5275, '3C':5341, '3D':5980,
                     "3' UTR":7264}

fermon_consensus = open('data/sequences/full_genome/fermon.fa').readlines()[1][1:]
mo_consensus = open('data/sequences/full_genome/18947.fa').readlines()[1][1:]

protein_old_new_codon = []
proteins = list(fermon_protein_pos.keys())
start_list = list(fermon_protein_pos.values())
print('Determining codon for each nucleotide')
protein_old_new_codon = []
proteins = list(fermon_protein_pos.keys())
start_list = list(fermon_protein_pos.values())
#for i in range(len(temp_melt)):
for i in range(10000):
    temp_row = melted_df.iloc[i]
    # Split the input mutation into original nucleotide, position, and new nucleotide
    genome_pos = temp_row['nucleotide position'] -1
    new_nuc = temp_row['variable']
    if temp_row['strain'] == 'Fermon':
        proteins = list(fermon_protein_pos.keys())
        start_list = list(fermon_protein_pos.values())
        consensus = fermon_consensus
    else:
        proteins = list(MO_protein_pos.keys())
        start_list = list(MO_protein_pos.values())
        consensus = mo_consensus
    nth_protein = sum(genome_pos > start_list) - 1
    protein_pos = genome_pos-start_list[nth_protein]-1
    #print(nth_protein, protein_pos)
    
    # Find the number of the gene amino acid
    aa_pos = (protein_pos) % 3 
    print(i, aa_pos, genome_pos)
    if nth_protein not in [0, 12]:
        # Return the amino acid number
        aa_num = ((protein_pos)//3)*3 + start_list[nth_protein]+1
        #print(i, aa_pos, nth_protein)

        # Find the old and new codon from this mutation
        
        old_codon = consensus[aa_num:aa_num+3]
        new_codon = old_codon[:aa_pos] + new_nuc + old_codon[aa_pos+1:]
        #print(old_codon, new_codon)
        # Find the old and new amino acid
        old_aa = str(Seq.Seq(old_codon).translate())
        new_aa = str(Seq.Seq(new_codon).translate())
        #print(old_aa, new_aa, old_codon, new_codon)
        # Separate the synonymous from non-synonymous mutation
        if old_aa != new_aa:
            non_synon = 'non-synonymous'
        else:
            non_synon = 'synonymous'
        protein_old_new_codon.append([consensus[genome_pos].upper(), proteins[nth_protein], 
                                      protein_pos, old_codon.upper(), new_codon.upper(), old_aa, new_aa, 
                                      non_synon, genome_pos, new_nuc])
    else:
        protein_old_new_codon.append([consensus[genome_pos].upper(), proteins[nth_protein], 
                                      protein_pos, 'NA', 'NA', 'NA', 'NA', 'NA', genome_pos, new_nuc])
protein_info_df = pd.DataFrame(protein_old_new_codon)    
protein_info_df.columns = ['WT nucleotide', 'protein', 'protein position', 'WT codon', 'mutant codon', 
                           'WT aa', 'mutant aa', 'non-synonymous', 'nucleotide position', 'variable']
protein_info_df = pd.DataFrame(protein_old_new_codon)    
protein_info_df.columns = ['WT nucleotide', 'protein', 'protein position', 'WT codon', 'mutant codon', 
                           'WT aa', 'mutant aa', 'non-synonymous', 'nucleotide position', 'variable']
print('Merging protein data with condition data')
full_table = pd.merge(left=temp_melt, right=protein_info_df, 
    left_on=['nucleotide position', 'variable'], right_on=['nucleotide position', 'variable'])
print('Writing to file')
full_table.to_csv('data/EV-D68_all_passage_data_v5.csv')
print('Done')
