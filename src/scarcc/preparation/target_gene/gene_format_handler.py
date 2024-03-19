"""This module contains functions to handle reading gene combinations from file and convert them to a standard format, able to handle checkerboard SG and DG"""

import ast
import pandas as pd
import itertools
import re

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

class GeneFormatHandler:
    def __init__(self, gene_str):
        self.is_Normal = 'Normal' in gene_str or '0.0' in gene_str
        self.gene_str = self.trim(gene_str)
        self.is_checkerboard_format = None
        self.DG = None
        self.first_gene = None
        self.second_gene = None

        self.SG = None
        self.XG = None
        self.classify_XG()
        self.set_fs_gene()

    @staticmethod
    def trim(gene_str):
        replacements = ["[A-Z]0_", "_(coculture|monoculture)"] # TODO: remove all id in species list
        return re.sub('|'.join(replacements),'', gene_str)

    def classify_XG(self):
        if '_' in self.gene_str:
            self.is_checkerboard_format = True
            genes, lvs = self.gene_str.split('_')
            self.first_gene, self.second_gene = [
                '_'.join(ele) for ele in zip(genes.split('.'), lvs.split('.'))] # g1_l1, g2_l2 
            if '0' in lvs:
                self.SG = self.gene_str
            else:
                self.is_DG = True
                self.DG = self.gene_str
        else:
            if '.' in self.gene_str:
                self.DG = self.gene_str
            else:
                self.SG = self.gene_str

    def set_fs_gene(self):
        if self.DG:
            if self.is_checkerboard_format:
                first_gene, second_gene = self.DG, self.DG
                first_gene = first_gene[:-1] + '0'
                second_gene = second_gene[:-3] + '0' + second_gene[-2:]
            else:
                first_gene, second_gene = self.DG.split('.')
            self.first_gene, self.second_gene = first_gene, second_gene
            self.SG = first_gene, second_gene

