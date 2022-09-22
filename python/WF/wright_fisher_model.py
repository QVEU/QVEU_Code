#!/usr/bin/env python3

# Benjamin Adam Catching
# NIH-NIAID-LVD
# Quantitative Virology and Evolution Unit

# Import modules
import numpy as np
import pandas as pd
#import scipy

import matplotlib.pyplot as plt
import seaborn as sns

# Define evolable class
class Agent:
    """
    An evolable object that represents one individual 'agent'

    attribute age: Number of cycles before this agent was created
    attribute fitness: How good is this agent at being selected
    attribute allele: Name of allele of this agent
    """

    def __init__(self, age, fitness, allele):
        self.age = age
        self.fitness = fitness
        self.allele = allele

def new_agents(pop_size, n_alleles, gen_age=0):
    """

    """

    # Define population size
    #pop_size = 100
    # Define dictionary of alleles
    num_allele = {}
    for i in range(n_alleles):
        num_allele[i] = chr(ord('@')+(i+1))


    # Initiate list for agents
    agents = []
    # Iterate and define each agent
    for i in range(pop_size):
        # Get the allele number
        temp_allele = i // (pop_size / n_alleles)
        agents.append(Agent(gen_age, 0.5, num_allele[temp_allele]))
    return agents

def rep_cycle(prev_agents, num_select):
    """
    Pick random agents from set of prior generation, no bias in selection

    param all_agents: List of agents to select from for new generation
    param num_select: How many agents to select

    return older_agents: New list of agents one cycle older
    """

    # Find the maximum number of iterations
    n_cycles = max(list(set([x.age for x in prev_agents])))
    # Get list of agents only with most recent iteration
    recent_agents = [x for x in prev_agents if x.age == n_cycles]

    """
    Potential issue with numpy.random.choice mutability, pick random manually
    """

    # Create set of agents, select and multiply for the agents
    new_agents = []
    # Get randomly selected list of agent alleles
    all_current_alleles = np.random.choice([x.allele for x in recent_agents],
                                            num_select, replace=True)
    # Select randomly from list
    for temp_allele in all_current_alleles:
        new_agents.append(Agent(n_cycles+1, 0.5, temp_allele))

    # Increment the age of the new all_agent
    older_agents = []
    counter = 0
    for a in new_agents:
        counter += 1
        a.age = n_cycles + 1
        older_agents.append(a)

    return older_agents

def replication_cycle(number_agents, number_alleles, number_cycles):
    """
    Neutral drift model of evolution, for a given number of alleles and cycles.

    param number_agents: number of agents per cycle to use
    """

    # Initiate population
    init_agents = new_agents(number_agents, number_alleles)
    all_agents = init_agents
    temp_iteration = init_agents

    # Set age counter
    n_cyc = 0
    # Iterate over the number of replication cycles
    while n_cyc != number_cycles:
        new_iteration = rep_cycle(temp_iteration, number_agents)
        all_agents += new_iteration
        temp_iteration = new_iteration

        n_cyc +=1

    # Format all agents into a dataframe
    run_df = pd.DataFrame([[x.age, x.allele] for x in all_agents])
    run_df.columns = ['age', 'allele']

    # Convert list of agents to DataFrame
    summary = []
    for a in set(run_df['allele']):
        for j in range(number_cycles+1):
            temp_df = run_df[(run_df['allele'] == a) &
            (run_df['age'] == j)]
            summary.append([j, a, len(temp_df)])
    # Summarize the number of each allele
    summary_df = pd.DataFrame(summary)
    summary_df.columns = ['age', 'allele', 'count']

    return run_df, summary_df


"""Test of Neutral Drift"""


# Define parameters
home_dir = '/Users/catchingba/Documents/Rosalind/'
n_agents = 10000
n_alleles = 4
n_cycles = 100

# Perform n_iterations replicates of evolution, save to dataframe
for rep in range(200):
    n_iterations = 0
    all_sum_df = pd.DataFrame()
    while n_iterations != 12:

        run_df, sum_df = replication_cycle(n_agents, n_alleles, n_cycles)
        num_columns = len(sum_df)
        sum_df['current run'] = n_iterations
        n_iterations += 1
        all_sum_df = pd.concat([all_sum_df, sum_df], ignore_index=False)
    out_file = home_dir + f'data/neutral/12reps_{rep}.csv'
    #print(rep)
    all_sum_df.to_csv(out_file)

# Plot example
sns.catplot(x='age', y='count', data=all_sum_df[all_sum_df['allele']=='A'],
            color='blue', hue='current run', kind='point', alpha=0.1)
plt.xlabel('passage number')
plt.ylabel('allele frequency')
plt.title('Allele Frequency over Time')
plt.tight_layout()
plt.savefig(home_dir+ '2022-07-21_neutral_allele_frequency.png', dpi=300)
#plt.show()
