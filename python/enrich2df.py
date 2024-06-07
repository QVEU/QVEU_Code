# Import packages
import numpy as np
import pandas as pd

# Open enrich2 data
enrich_df = pd.read_csv('../../data/sequences/enrich_3D_DMS/enrich_output/tsv/WT_and_ribavirin_exp/main_identifiers_scores_shared_full.tsv',  delimiter='\t')

"""RUN ONLY ONCE for the output of Enrich2"""

# Define the number of conditions run in this analysis
num_conditions = 4
#
conditions = [x.split('.')[0] for x in enrich_df.columns.to_numpy()]
replicates = enrich_df.iloc[0].to_numpy()
values = enrich_df.iloc[1].to_numpy()
new_columns = ['_'.join(np.array([conditions, replicates, values])[:, i]) for i in range(num_conditions* 6 + 1)]

enrich_df = enrich_df.iloc[4:]
enrich_df.columns = ['value'] + new_columns[1:]

enrich_df['mutant'] = [x.split('(')[1][:-1] for x in enrich_df['value']]
enrich_df['aa position'] = [int(x[1:-1]) for x in enrich_df['mutant']]
enrich_df['aa'] = [x[-1] for x in enrich_df['mutant']]

# Output file
enrich_df.to_csv('../../data/sequences/enrich_3D_DMS/enrich_df.csv')