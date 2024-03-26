"""Checkerboard Alpha Finder"""

import numpy as np
import pandas as pd
import itertools
import os
import ast

from .monoculture import MonocultureAlphaFinder

def read_normalized_growth(normalized_growth_file_path, model_list):
    """
    Read normalized growth file that is readily usable by the class CheckerboardAlphaFinder
    Requirement of the input file, each cell stores a list of normalized biomass values for a given model and gene.

    The function performs the following steps:
    1. Reads the CSV file at the provided path into a pandas DataFrame.
    2. Replaces the column names with the corresponding model objects from the provided model list.
    3. Converts string representations of lists in the DataFrame into actual list objects.

    Parameters
    ----------
    normalized_growth_file_path (str): Path to the CSV file containing normalized growth data.
    model_list (list): List of model objects, used to match model IDs in the CSV file.

    Returns
    -------
    pd.DataFrame: A DataFrame containing the normalized growth data, with model objects as column names 
        and list od normalized biomass as cell values.
    """
    def sub_id_to_model(model_id):
        for model in model_list:
            if model.id == model_id:
                return model
        return model_id

    def convert_to_list(l):
        return l if isinstance(l, list) else ast.literal_eval(l)

    normalized_growth = pd.read_csv(normalized_growth_file_path, index_col=0)
    normalized_growth.columns = [sub_id_to_model(id)
                                    for id in normalized_growth.columns]
    normalized_growth = normalized_growth.applymap(convert_to_list)
    return normalized_growth

class CheckerboardAlphaFinder():
    """
    Checkerboard Alpha Finder that calls MonocultureAlphaFinder to find the alpha values for each gene in each species

    Parameters
    ----------
    response_record: dict, optional
        A dictionary containing the normalized growth data for each gene in each species.
    normalized_growth: pd.DataFrame
        Output from read_normalized_growth function, containing the normalized target biomass to perform search on.
    data_directory: str
        Path of Data directory to read normalized growth data and save the alpha table.
    **maf_kwargs: dict
        Additional keyword arguments to pass to the MonocultureAlphaFinder class.

    Attributes
    ----------
    response_record (dict): A dictionary containing the normalized growth data for each gene in each species.
        The dictionary should have the following structure, for gene i and species j:
        {
            gene_i: {
                species_j: {
                    nbiomass_x: [list of normalized biomass values],
                    response: {normalized_biomass1: alpha1, normalized_biomass2: alpha2, ...}
                },
                ...
            },
            ...
        }
        Each 'response' dictionary is with keys normalized biomass from 0 to 1 rounded to nearest 0.05, with values as alpha that gives the normalized biomass.
        This nested dictionary is gradually updated in MonocultureAlphaFinder, keeping all searched results.
    target_gene_list (list): List of gene names in the response_record.
    species_list (list): List of cobra.model in the response_record.
    n_levels (int): Number of normalized biomass levels. ! Assumed all target genes were simulated with same number of levels
    checker_lvs (pd.DataFrame): DataFrame containing all possible pairs of normalized biomass levels.
    result_dict (dict): Dictionary containing the alpha table rows for each lv_pairs.
    data_directory (str): Path of Data directory to read normalized growth data and save the alpha table.
    maf_kwargs (dict): Additional keyword arguments to pass to the MonocultureAlphaFinder class.
    alpha_table (pd.DataFrame): DataFrame containing the alpha values for each gene in each species.
        The DataFrame has the following columns:
        - Gene_inhibition: Name of the single target gene.
        - lv_pairs: Tuple of two ordinal numbers representing the concentration levels of each drug
        - E0_alpha: Alpha value for E. coli.
        - E0_lv: Ordinal number of the first concentration level for E.coli.
        - E0_normalized_biomass: Normalized biomass value for E.coli under perturbation parameters alpha
        - S0_alpha: Alpha value for S. enterica
        - S0_lv: Ordinal number of the first concentration level for S. enterica.
        - S0_normalized_biomass: Normalized biomass value for S. enterica under perturbation parameters alpha
    """

    def __init__(self, response_record=None, normalized_growth=None, data_directory = None, **maf_kwargs) -> None:
        self.response_record = self.set_response_record_format(response_record, normalized_growth)
        self.target_gene_list = list(self.response_record.keys())
        self.species_list = list(self.response_record[self.target_gene_list[0]].keys())
        
        # checkerboard levels
        nbiomass_x = self.response_record[self.target_gene_list[0]][self.species_list[0]]['nbiomass_x']
        self.n_levels = len(nbiomass_x)
        self.checker_lvs = self.construct_checker_lvs()
        self.result_dict = {}
        self.data_directory = data_directory
        # default kwargs passed to MonocultureAlphaFinder
        self.maf_kwargs = maf_kwargs
        self.set_default_kwargs('acceptance_threshold_upper', .995)
        self.set_default_kwargs('acceptance_threshold_lower', 1.001)
        self.set_default_kwargs('precision', 3)
        # output
        self.alpha_table = None

    @staticmethod
    def set_response_record_format(response_record, normalized_growth):
        """
        construct response_record to initialize MonocultureAlphaFinder
        
        Parameters
        ----------
        response_record: dict, optional
            A dictionary containing the normalized growth data for each gene in each species.
        normalized_growth: pd.DataFrame
            Output from read_normalized_growth function, containing the normalized target biomass to perform search on.

        Returns
        -------
        response_record: dict

        """
        if response_record is None:
            response_record = normalized_growth.to_dict(orient='index')
            response_record = {
                current_gene: {
                    model: {'nbiomass_x': nbiomass_list}
                    for model, nbiomass_list in model_dict.items()
                }
                for current_gene, model_dict in response_record.items()
            }
        return response_record

    @staticmethod
    def get_new_response_record(**kwargs):
        """get MonocultureAlphaFinder object to update response_record
        Parameters
        ----------
        **kwargs: dict
            Keyword arguments to pass to the MonocultureAlphaFinder class.
            
        Returns
        -------
        MonocultureAlphaFinder: MonocultureAlphaFinder object
        """
        maf = MonocultureAlphaFinder(**kwargs)
        maf.find_feasible_alpha()
        return maf

    def set_default_kwargs(self, key, default_value):
        """Set default value for maf_kwargs if key not in maf_kwargs"""
        if key not in self.maf_kwargs:
            self.maf_kwargs[key] = default_value

    def fill_dict(self, model, current_gene, nbiomass_x):
        """response_record get updated with new desired normalized biomass(nbiomass_x) that needs to be searched"""
        nbiomass_lv = float(format(nbiomass_x, '.2f'))
        inner_dict = self.response_record[current_gene][model].setdefault('response', {})
        visited= inner_dict.keys()
        if float(format(nbiomass_lv, '.2f')) not in visited:
            if inner_dict.get(nbiomass_lv):
                closest_val = min(visited, key=lambda x: abs(x - nbiomass_x))
                search_alpha = self.response_record[current_gene][model]['response'][closest_val]
                search_alpha = search_alpha*1.1 if closest_val > nbiomass_x else search_alpha*0.9

            else:
                search_alpha = 1.2
            maf = self.get_new_response_record(
                model=model, search_alpha=search_alpha,
                current_gene=current_gene, target_normalized_biomass=nbiomass_lv,
                response_record=self.response_record, **self.maf_kwargs)
            return maf
        return self.response_record

    def construct_checker_lvs(self):
        """Construct all possible pairs of normalized biomass levels
        
        Return:
        -------
        pd.DataFrame: DataFrame containing all possible pairs of normalized biomass levels and level for each species
            will be used in constructing the alpha table.
        """
        checker_lvs = pd.DataFrame()
        MIC_levels = np.arange(0, self.n_levels)

        checker_lvs['lv_pairs'] = list(itertools.product(MIC_levels, MIC_levels))
        checker_lvs['lv1'] = checker_lvs['lv_pairs'].apply(lambda x: x[0])
        checker_lvs['lv2'] = checker_lvs['lv_pairs'].apply(lambda x: x[1])
        return checker_lvs

    def match_alpha_lvs(self, current_gene):
        """Pull the correct alpha values for each species in the response_record"""
        species_alpha_table = {}
        for Species in self.species_list:
            record_in_gene_species = self.response_record[current_gene][Species]
            temp_dict = {
                ith: {
                    'alpha': record_in_gene_species['response'][float(format(nbiomass, '.2f'))]['search_alpha'],
                    'lv': ith,
                    'Gene_inhibition': current_gene,
                    'normalized_biomass': nbiomass}
                        for ith, nbiomass in enumerate(record_in_gene_species['nbiomass_x'])}
            species_alpha_table[Species] = pd.DataFrame.from_dict(temp_dict, orient='index').set_index('Gene_inhibition').add_prefix(f'{Species.id}_')
        return species_alpha_table
    
    def process_record_in_gene(self, current_gene):
        """Merge the alpha values for each species in the response_record to get the row for each lv_pairs"""
        outer_dict = self.match_alpha_lvs(current_gene)
        for sp1, sp2 in itertools.combinations(outer_dict.keys(), 2):
            self.result_dict[current_gene, sp1, sp2] = (
                self.checker_lvs.merge(outer_dict[sp1], left_on='lv1', right_on=f'{sp1.id}_lv')
                .merge(outer_dict[sp2], left_on='lv2', right_on=f'{sp2.id}_lv', suffixes=(f'_{sp1}', f'_{sp2}'))
                .drop(columns=['lv1', 'lv2']))
            self.result_dict[current_gene, sp1, sp2]['Gene_inhibition'] = current_gene
            self.result_dict[current_gene, sp1, sp2].set_index('Gene_inhibition', inplace=True)

    def get_checkerboard_alpha_table(self):
        """Get the final alpha values ready for checkerboard simulation"""
        for current_gene, sub_dict in self.response_record.items():
            for model, ss_dict in sub_dict.items():
                for nbiomass in sorted(ss_dict['nbiomass_x']):
                    maf = self.fill_dict(model, current_gene, nbiomass)
                    self.response_record = maf.response_record  if isinstance(maf, MonocultureAlphaFinder) else maf
                    alpha_lvs = ss_dict.setdefault('alpha_lvs', {})
                    IC_lv = float(format(nbiomass, '.2f'))
                    alpha_lvs.update({nbiomass: ss_dict['response'][IC_lv]['search_alpha']})

        _ = [self.process_record_in_gene(current_gene) for current_gene in self.target_gene_list]
        self.alpha_table = pd.concat(self.result_dict.values(), axis=0).sort_values('lv_pairs')

        if self.data_directory is not None:
            self.alpha_table.to_csv(os.path.join(self.data_directory, 'alpha_table_checkerboard.csv'))
            print(f'Checkerboard alpha table saved to {self.data_directory}.')
        return self.alpha_table