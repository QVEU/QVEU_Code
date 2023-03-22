# Benjamin Adam Catching
# 2023-03-02
# NIH-NIAID-LVD
# Quantitative Virology and Evolution Unit

"""
Oligocottus maculosus (Tidepool sculpin)
From aligned, phred-score = 20 cutoff, Q20.txt files, call each consensus
mutation in each replicate of a 8-replicate, 12-condition passaging experiment
"""

# Import modules
import numpy as np
import pandas as pd
import glob

import matplotlib.pyplot as plt
import seaborn as sns

from Bio import Seq

"""====Define functions===="""

# Define function for extracting consensus and >50% mutants
def Q20_to_consensus(q20_df, count_thresh=20):
    """
    Read through dataframe of position nucleotide counts
    
    param: q20_df, Pandas DataFrame encompassing nucleotide counts

    param: count_thresh, threshold of coverage at each position required
    to declare the consensus

    return: consensus_list, list of consensus nucleotides
    
    return: sums, list of total number of nucleotides at each position
    """
    
    # Dictonary to convert position to nucleotide
    pos_to_nuc = {0:'A', 1:'C', 2:'G', 3:'T'}
    # Enumerate list of max nucleotides 
    consensus_list = []
    sums = []
    for i in range(len(q20_df)):
        if sum(q20_df.iloc[i]) >= count_thresh:
            new_nuc = pos_to_nuc[list(q20_df.iloc[i]).index(max(list(q20_df.iloc[i])))]
        else:
            new_nuc = 'N'
        sums.append(q20_df.iloc[i].sum())
    return consensus_list, sums


def coverage_consensus_96(Q20_files):
	"""
    Iterate over all Q20 input files to return a consensus sequence for each 
    Q20 file

    param: Q20_files, a list of Q20 text files, foreach  to be each found a 
    consensus sequence

    return: consensus_list, list of consensus sequences
    
    return: sums, list of lists of total number of nucleotides at each position
	"""

	consensus_list = []
	coverage = []

	# Iterate over the 
	for i in range(96):
        i += 1
        if i % 8 == 0:
        # Define next file
        temp_file = Q20_files[i]
        # Open file
        temp_df = pd.read_csv(temp_file, delimiter='\t')
        temp_df.columns = ['A', 'C', 'G', 'T', 'N']
    
        # Get consensus
        temp_consensus = Q20toconsensus(temp_df)
        consensus_list.append(''.join(temp_consensus[0]))
        coverage.append(temp_consensus[1])

    return consensus_list, coverage


def mut_count(list_of_list_of_muts):
    """
    Count the number of occurances of a mutation and store in a dictionary

    param: list_of_list_of_muts, a list containing a list of mutations

    return: mut_count, a dictionary of each mutation and the number of times
    it was observed
    """
    mut_count = {}
    for rep in list_of_list_of_muts:
        for mut in rep:
            if mut in mut_count.keys():
                mut_count[mut] += 1
            else:
                mut_count[mut] = 1
    return mut_count


def consensus_mutations(offset=False, consensus_list):
    """
    From a list of 96 consensus sequences call the mutants when compared to a 
    reference sequence

    param: offset, boolean value of whether to start with 0 or 8 as starting
    position of the batch (96 samples batched as 12 batches of 8)

    param: consensus_list, list of strings representing all 96 called consensus
    sequences

    return: ref_mut_df, DataFrame of the conditions each mutation with the 
    number of observations per occurance
    """

    # Define the list of consensus 
    batches = []
    # Offset the count of the barcode number by groups of 8
    if offset == True:
        batch_offset = 8
        ref_length = Fermon_length
        ref_consensus = Fermon_consensus
        strain = 'Fermon'
    else:
        batch_offset = 0
        ref_length = MO_length
        ref_consensus = MO_consensus
        strain = 'US/MO/14-18947'
    for i in range(6):
        temp_batch_muts = []

        for j in range(8):
            
            j += batch_offset + (i * 16)   
            temp_consensus = consensus_list[j]
            # Store since lineage mutations
            temp_muts = []
            for k in range(ref_length):
                temp_nuc = temp_consensus[k]
                ref_nuc = ref_consensus[k]
                if temp_nuc != ref_nuc and 'N' not in [temp_nuc, ref_nuc]:
                    temp_muts.append(''.join([ref_nuc, str(k+1), temp_nuc]))
            temp_batch_muts.append(temp_muts)
        batches.append(temp_batch_muts)
        
    # Sort and display US/MO/2014-18947
    batch_sets = [mut_count(x) for x in batches]

    ref_conditions = [strain+', SH-SY5Y, 33C', strain+', SH-SY5Y, 37C', 
                      strain+', RD, 33C', strain+', RD, 37C', 
                      strain+', A549, 33C', strain+', A549, 37C']
    all_mutations = []
    for i in range(6):
        condition = ref_conditions[i].split(', ')
        print(condition)
        for j in batch_sets[i].items():
            all_mutations.append(condition + [j[0], j[1]])

    ref_mut_df = pd.DataFrame(all_mutations)
    ref_mut_df.columns = ['strain', 'cell', 'temperature', 
                          'mutation', 'number']

    return ref_mut_df


def mutation_info(gene_dict, temp_mutations, ref_cons):
    """
    From a dictionary of gene positions, a list of mutations, and the 
    reference consensus genome,return a DataFrame of the information of 
    each mutation
    
    param: gene_dict, a dictionary of each element in the genome and their 
    start position
    
    param: temp_mutations, a list of mutations in the form X###X

    param: ref_cons, the reference consensus sequence

    return: mut_df, a DataFrame of mutations with their nucleotide position, 
    protein, aa mutation, old codon, new codon, nucleotide mutation, and 
    whether the mutation is synonymous or not
    """
    
    
    # Array of elements, including 5' UTR to account for 0-th counting
    gene_list = np.array(["5' UTR"] + list(gene_dict.keys()))
    # Array of start positions, to quickly determine which protein
    start_list = np.array(list(gene_dict.values()))
    # List of each gene
    sequences = [ref_cons[gene_dict[gene_list[i+1]]:gene_dict[gene_list[i+2]]] for i in range(11)]
    
    # List to store mutation infomation
    high_frequency_muts = []
    # The number of mutations to iterate through
    num_mutations = len(temp_mutations)
    
    # Iterate over input mutations
    for i in range(num_mutations):
        # Define the mutation to be examined
        mut_of_interest = temp_mutations[i]
        # Split the input mutation into original nucleotide, position, and new nucleotide
        temp_mut = [mut_of_interest[0], int(mut_of_interest[1:-1]), mut_of_interest[-1]]
        # Quickly find the element the mutation is in
        nth_protein = sum(temp_mut[1] > start_list)
        
        # Define the genome and gene position
        genome_pos = temp_mut[1]
        protein_pos = temp_mut[1]-start_list[nth_protein-1]
        # Find the number of the gene amino acid
        aa_pos = (protein_pos-1) % 3

        if nth_protein not in [0, 12]:
            # Return the amino acid number
            aa_num = (protein_pos+2) //3 
            
            # Find the old and new codon from this mutation
            old_codon = sequences[nth_protein-1][aa_num*3-3:aa_num*3]
            new_codon = old_codon[:aa_pos] +mut_of_interest[-1] + old_codon[aa_pos+1:]
            
            # Find the old and new amino acid
            old_aa = str(Seq.Seq(old_codon).translate())
            new_aa = str(Seq.Seq(new_codon).translate())
            # Separate the synonymous from non-synonymous mutation
            if old_aa != new_aa:
                non_synon = 'non-synonymous'
            else:
                non_synon = 'synonymous'
            high_frequency_muts.append([int(temp_mut[1]), 
                gene_list[nth_protein], old_aa + str(aa_num) + new_aa, 
                old_codon, new_codon, mut_of_interest, non_synon])
                
    mut_df = pd.DataFrame(high_frequency_muts)
    mut_df.columns = ['nucleotide position', 'protein', 'aa mutation', 
    'old codon', 'new codon', 'nucleotide mutation', 'synonymous or not']

    return mut_df.sort_values(by='nucleotide position', ascending=True)


"""
Creation of surfseq files
"""


# Extract locations of all Q20 files
global_dir = '../../data/sequences/'
P1_files = sorted(glob.glob(global_dir + '2023-02-06_passage_1/Q20/*.txt'))
P2_files = sorted(glob.glob(global_dir + '2023-02-27_passage_2/Q20/*.txt'))
P3_files = sorted(glob.glob(global_dir + '2023-03-14_passage_3/Q20/*.txt'))

# Find locations of the sequencing from the initial virus (US/MO/2014-18947=MO, 
# Fermon=Fermon)
base_dir = global_dir+'2023-01-20_ARTIC/mk1b/ultra_slow_basecalling/aligned/'
MO_P0_Q20_loc = base_dir + '18947_aligned_Q20.txt'
Fermon_P0_Q20_loc = base_dir + 'fermon_aligned_Q20.txt'

# Import data as dataframe
MO_P0_Q20_df = pd.read_csv(MO_P0_Q20_loc, delimiter='\t')
MO_P0_Q20_df.columns = ['A', 'C', 'G', 'T', 'N']
Fermon_P0_Q20_df = pd.read_csv(Fermon_P0_Q20_loc, delimiter='\t')
Fermon_P0_Q20_df.columns = ['A', 'C', 'G', 'T', 'N']

# Define reference lengths
MO_length = 7370
Fermon_length = 7364
MO_consensus = ''.join(Q20toconsensus(MO_P0_Q20_df)[0])[:MO_length]
Fermon_consensus = ''.join(Q20toconsensus(Fermon_P0_Q20_df)[0])[:Fermon_length]

"""Passage 3 data"""
P3_consensus_list = coverage_consensus_96(P3_files)
P3_MO_df = consensus_mutations(offset=False, P3_consensus_list[0])
P3_fermon_df = consensus_mutations(offset=True, P3_consensus_list)

MO_muts = mutation_info(MO_protein_pos, list(set(P3_MO_df['mutation'])), 
    MO_consensus)
P3_MO_df = P3_MO_df.merge(MO_muts, left_on='mutation', 
    right_on='nucleotide mutation').sort_values(by='nucleotide position')



"""Plotting data"""

"""
# Plotting Coverage
sns.set_context('talk')

conditions = ['US/MO/14-18947\nSH-SY5Y, 33C', 'Fermon\nSH-SY5Y, 33C', 
              'US/MO/14-18947\nSH-SY5Y, 37C', 'Fermon\nSH-SY5Y, 37C',
              'US/MO/14-18947\nRD, 33C', 'Fermon\nRD, 33C', 
              'US/MO/14-18947\nRD, 37C', 'Fermon\nRD, 37C',
             'US/MO/14-18947\nA549, 33C', 'Fermon\nA549, 33C', 
             'US/MO/14-18947\nA549, 37C', 'Fermon\nA549, 37C']
P3_coverage = P3_consensus_list[1]
with sns.axes_style('whitegrid'):
    fig, ax = plt.subplots(8, 12, figsize=(36, 24), sharey=True, sharex=True)
    for i in range(12):
        # i is the condition
        for j in range(8):
            # j is the replicate
            
            # Store current Q20 list
            temp_Q20 = P3_coverage[i*8+j]
            ax[j, i].plot(np.linspace(1,len(temp_Q20), len(temp_Q20)), temp_Q20, 
            color='darkslategray')
            ax[j, i].set_yscale('log')
            ax[j, i].set_xticks([0,2500,5000,7500])
            ax[j, 0].set_ylabel('coverage')
        if j == 7:
            ax[j, i].set_xlabel('genome\nposition')
            ax[7, i].tick_params(axis='x', rotation=45)
        ax[0, i].set_title(conditions[i])
    plt.savefig('../../2023-03-15_P3_coverage.png', dpi=300)
"""


"""
# Showing consensus mutations
# Divide the number of subplots into 3 by 2
fig, ax = plt.subplots(2, 3, figsize=(20, 10), sharey=True, sharex=True)

cell_types = list(set(P3_MO_df['cell']))
temps = list(set(P3_MO_df['temperature']))

for j in range(3):
    for k in range(2):
        temp_df = P3_MO_df[(P3_MO_df['cell'] == cell_types[j]) & \
        (P3_MO_df['temperature'] == temps[k])]
        temp_count = list(temp_df['number'])
        temp_pos = list(temp_df['nucleotide position'])
        syn_nonsyn = [int(x == 'synonymous') for x in list(temp_df['synonymous or not'])]
        protein = list(temp_df['protein'])
        temp_mut = list(temp_df['aa mutation'])
        colors = ['navy', 'firebrick']
        #print(syn_nonsyn)
        for i in range(len(temp_count)):
            ax[k, j].plot(temp_pos[i], temp_count[i], 'o', color=colors[syn_nonsyn[i]], alpha=0.25)
ax[0,0].set_ylabel('')
"""