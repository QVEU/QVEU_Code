# Benjamin Adam Catching
# 2022-10-24
# NIH-NIAID-LVD-QVEU
# EV-D68 Thermal Adaptation Project

"""
Read in lines of the SAM file and parse each aligned string into the complete
sequence with spaces at gaps
"""

# Import packages
import sys
import numpy as np

# Define dictionaries
phred_word = '!"'+"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
letter2phred = {}
for i, p in enumerate(phred_word):
    letter2phred[p] = i

seq2pos = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}

# Read in sequence
# This needs to be versitile for different reads
WT_file = open(str(sys.argv[1]))
WT_lines = WT_file.readlines()
WT_file.close()
WT_seq = ''.join(WT_lines[1:])
WT_len = len(WT_seq)
#print(WT_seq, WT_len)

# Define the directory
temp_file = str(sys.argv[2])

# Define where to store the nucleotide occurances
q8_file = np.zeros((WT_len, 5))

# Open file
file = open(temp_file, 'r')

# Define the number of lines in the file
num_lines = sum(1 for _ in file)
print(num_lines)
oneper = int(num_lines /100)
file.close()
file = open(temp_file, 'r')

# Read a chunk of data
file.readline()
file.readline()

# Initialize file for nucleotide counts
count_file = open(str(sys.argv[3]), 'w')

# Initialize the number of reads
num_reads = 0
# Initialize the number of matched reads
virus_reads = 0
for k in range(num_lines-2):

    num_reads += 1
    # Parse out data from line
    temp_line = file.readline()
    #print(temp_line)
    # Separate the different data from the line
    temp_data = temp_line.split('\t')
    #print(temp_data)
    # Define sequence
    temp_start = temp_data[3]
    temp_sequence = temp_data[9]
    temp_CIGAR = temp_data[5]
    temp_phred = temp_data[10]

    """
    Reconstruct the sequence
    """
    if temp_start != 0:
        virus_reads += 1

        # Define the new sequence string
        new_aligned_seq = 'N' * (int(temp_start) - 1)
        new_aligned_fred = '!' * (int(temp_start) - 1)

        # Split CIGAR into each section
        split_CIGAR = [x for x in temp_CIGAR]
        # Define the prior number of nucleotides by list of numbers
        temp_number = []
        # Define the current working position
        current_pos = 0
        # Iterate over the integer and value of the CIGAR string`
        for i, x in enumerate(split_CIGAR):
            # If the value is a letter, intepret the preceeding integer
            if x in ['I', 'S', 'M', 'D', 'H']:
                integer_pair = int(''.join(temp_number))
                #print('  '+x)

                # If match add the previous values
                if x == 'M':
                    end_pos = current_pos+integer_pair
                    new_aligned_seq += temp_sequence[current_pos:end_pos]
                    new_aligned_fred += temp_phred[current_pos:end_pos]
                    # If the letter corresponds to a mismatch, add gap of spaced
                elif x == 'H':
                    # Write gap of spaces
                    temp_gap = 'N' * integer_pair
                    temp_score = '!' * integer_pair
                    # Add the gap to the string
                    new_aligned_seq += temp_gap
                    new_aligned_fred += temp_score
                elif x == 'D':
                    temp_gap = 'N' * integer_pair
                    temp_score = '!' * integer_pair
                    integer_pair =0
                    new_aligned_seq += temp_gap
                    new_aligned_fred += temp_score


                # Update integer list and current position
                temp_number = []
                # Add the gap to the string
                current_pos += integer_pair
            else:
                #print(x)
                temp_number.append(x)
        for j in range(len(new_aligned_seq)-1):
            temp_nt = new_aligned_seq[j]

            #print(j, new_aligned_fred[j], letter2phred[new_aligned_fred[j]])
            temp_phred = letter2phred[new_aligned_fred[j]]
            if j < WT_len and temp_phred >= 20:
                temp_nt = new_aligned_seq[j]
                q8_file[j, seq2pos[temp_nt]] += 1

file.close()
for line in q8_file:
    count_file.write('\t'.join([str(int(x)) for x in line]) + '\n')
count_file.close()
print('total reads: ' + str(num_reads) + '\n' + 'viral reads: ' + str(virus_reads))
#print(q8_file)
