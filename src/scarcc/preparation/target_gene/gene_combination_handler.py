"""This module contains functions to handle reading gene combinations from file and convert them to a standard format, able to handle checkerboard SG and DG"""

import ast
import pandas as pd
import itertools

def get_DG_list(filepath, n_combos=None):
    DG_list = pd.read_csv(filepath, header=None)
    DG_list = [ast.literal_eval(i) for i in DG_list.iloc[:,0].to_list()]
    if n_combos is not None:
        DG_list = DG_list[:n_combos]
    if 'Gene_inhibition' in DG_list:
        DG_list.remove('Gene_inhibition')
    return DG_list

def get_SG_list(DG_list):
    if isinstance(DG_list[0], str) and '.' in DG_list[0]:
        DG_list = [ele.split('.') for ele in DG_list]
    SG_list = set(itertools.chain(*DG_list))
    return list(SG_list)

def generate_all_combinations(SG_list):
    if isinstance(SG_list, str):
        raise ValueError('SG_list must be a list of target genes')
    potential_combos = list(itertools.combinations(SG_list, 2))
    print(f'There are {len(potential_combos)} potential gene combinations')
    return potential_combos