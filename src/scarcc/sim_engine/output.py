"""Module for producing output data from simulation and saving in data directory"""

import pandas as pd
import os
import itertools
import logging

from scarcc.data_analysis import convert_po_col
from scarcc.data_analysis.growth.growth_rate import get_growth_rate_df
from scarcc.preparation.target_gene.gene_format_handler import GeneFormatHandler

logger = logging.getLogger(__name__)

def merge_single_gene_gr(sg_df, dg_df):
    """Merge single gene growth rate data from either first or second gene in drug combination as supplied in dg_df index"""
    dg_df = pd.concat([sg_df.loc[['Normal']], dg_df]).rename(index={'Normal': 'Normal.Normal'})
    dg_df['First_gene'], dg_df['Second_gene'] = list(zip(*dg_df.index.map(lambda x: GeneFormatHandler(x).SG)))
    dg_df.rename(index={'Normal.Normal': 'Normal'}, inplace=True)
    dg_df = (dg_df
                .merge(sg_df, left_on='First_gene', right_index=True, suffixes=['','_First_gene'])
                .merge(sg_df, left_on='Second_gene', right_index=True, suffixes=['','_Second_gene'])
                )
    return dg_df

def add_additive_and_drug_comb_response(df, additive_threshold=0.05): # after normalization of DG_gr
    """Add columns for classification of drug combination response based on sign of difference between observed and predicted values.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing growth rate data for drug combination experiments.
    additive_threshold : float, optional
        Threshold for classifying drug combination response as additive. Default is 0.05.
        
    Returns
    -------
    pd.DataFrame
        DataFrame with added columns for classification of drug combination response.
    """
    for sp_cul in itertools.product(['E0', 'S0'], ['coculture', 'monoculture']):
        sp_cul = '_'.join(sp_cul)
        if f'{sp_cul}_First_gene' in df.columns:
            df[f'Predicted_additive_effect_{sp_cul}'] = df[f'{sp_cul}_First_gene'] * df[f'{sp_cul}_Second_gene']
            # df[f'op_diff_{sp_cul}'] = df[sp_cul] - df[f'Predicted_additive_effect_{sp_cul}'] # TODO: call convert_op_col
            df[f'po_diff_{sp_cul}'] = df[f'Predicted_additive_effect_{sp_cul}'] - df[sp_cul] # predicted - observed
            df[f'Drug_comb_effect_{sp_cul}'] = convert_po_col(df[f'po_diff_{sp_cul}'], additive_threshold=additive_threshold)
    return df

def reorder_columns(df):
    """Reorder columns in DataFrame if they are present"""
    column_order = ['E0_coculture', 'S0_coculture',
    'Predicted_additive_effect_E0_coculture',
    'Predicted_additive_effect_S0_coculture', 'E0_monoculture',
    'S0_monoculture', 'Predicted_additive_effect_E0_monoculture',
    'Predicted_additive_effect_S0_monoculture', 'First_gene',
    'E0_monoculture_First_gene', 'S0_monoculture_First_gene', 'Second_gene',
    'E0_coculture_Second_gene', 'S0_coculture_Second_gene',
    'E0_monoculture_Second_gene', 'S0_monoculture_Second_gene',
    'po_diff_E0_coculture', 'po_diff_S0_coculture',
    'po_diff_E0_monoculture', 'po_diff_S0_monoculture', 
    'Drug_comb_effect_E0_coculture', 'Drug_comb_effect_S0_coculture',
    'Drug_comb_effect_E0_monoculture', 'Drug_comb_effect_S0_monoculture']
    
    return df.reindex(columns=[col for col in column_order if col in df.columns])

class DrugCombinationEffectClassification:
    """Class for adding columns for classification of drug combination response with corresponding growth rate data in single target gene given
    
    Parameters
    ----------
    dfs_in_SG_layer : dict
        Dictionary containing data frames for single gene experiments.
    dfs_in_DG_layer : dict
        Dictionary containing data frames for drug combination experiments.
    additive_threshold : float, optional
        Threshold for classifying drug combination response as additive. Default is 0.05.

    Attributes
    ----------
    dfs_in_SG_layer : dict
        Dictionary containing data frames for single gene experiments.
    dfs_in_DG_layer : dict
        Dictionary containing data frames for drug combination experiments, if not None, add drug combination response classification.
    pure_dg_gr : pd.DataFrame
        DataFrame containing growth rate data for drug combination experiments before merging with single gene growth rate data.
    additive_threshold : float
        Threshold for classifying drug combination response as additive.
    """
    def __init__(self, dfs_in_SG_layer, dfs_in_DG_layer, additive_threshold=0.05):
        self.dfs_in_SG_layer = dfs_in_SG_layer
        self.dfs_in_DG_layer = dfs_in_DG_layer
        self.pure_dg_gr = None
        self.additive_threshold = additive_threshold

    def get_gr_df(self):
        """Add gr_df to the container dict from biomass

        1. identify Normal case as normalizing factor
        2. merge single gene gr data to dg gr data
        3. reorder columns
        """
        sg_df = get_growth_rate_df(self.dfs_in_SG_layer['biomass'])
        if 'Normal' not in sg_df.index:
            Normal_row = sg_df.loc[[index for index in sg_df.index if '0.0' in index][0]]
            sg_df.loc['Normal'] = Normal_row
        self.dfs_in_SG_layer['growth_rate'] = reorder_columns(sg_df)
        if self.dfs_in_DG_layer is not None:
            dg_df = get_growth_rate_df(self.dfs_in_DG_layer['biomass'])
            self.pure_dg_gr = dg_df.copy() # This is before merge, because joining sg_df in non-normalized and normalized

            dg_df = merge_single_gene_gr(sg_df, dg_df)
            self.dfs_in_DG_layer['growth_rate'] = reorder_columns(dg_df)

    def get_normalized_gr_df(self):
        """Add normalized_gr_df to the container dict from growth_rate

        1. normalize growth rate data by dividing by Normal row
        2. merge single gene gr data to dg gr data
        3. reorder columns
        """
        sg_df = self.dfs_in_SG_layer['growth_rate']
        Normal_row = sg_df.loc['Normal']
        sg_df = sg_df.div(Normal_row)
        self.dfs_in_SG_layer['normalized_growth_rate'] = reorder_columns(sg_df)
            
        if self.dfs_in_DG_layer is not None:
            dg_df = self.pure_dg_gr
            dg_df = dg_df.div(Normal_row)
            dg_df = merge_single_gene_gr(sg_df, dg_df) # merged one is normalized
            dg_df = add_additive_and_drug_comb_response(dg_df, additive_threshold=self.additive_threshold)
            self.dfs_in_DG_layer['normalized_growth_rate'] = reorder_columns(dg_df)
            self.dfs_in_DG_layer['drug_response_classification'] = dg_df.filter(regex='Drug_comb_effect_')

class MethodDataFiller:
    """Class for filling df_container dict with tuple (method, XG) as keys

    Parameters
    ----------
    df_container : dict
        Dictionary containing data frames for single gene and drug combination experiments.
    data_directory : str
        Path to directory for saving output data.
    additive_threshold : float, optional
        Threshold for classifying drug combination response as additive. Default is 0.05.

    Attributes
    ----------
    df_container : dict
        Dictionary containing data frames for single gene and drug combination experiments.
    data_directory : str
        Path to directory for saving output data.
    methods : set
        Set of methods in df_container.
    XGs : set
        Set of XGs in df_container.
    additive_threshold : float
        Threshold for classifying drug combination response as additive.
    """
    # TODO: handler for not supplying SG biomass, OR handle for the first construction of container dict
    # TODO: if DG none then no cal any for DG OR ensure both SG and DG exist
    def __init__(self, df_container, data_directory, additive_threshold=0.05): #? make a base accept only sg, dg
        self.df_container = df_container
        self.data_directory = data_directory
        self.methods = set([key[0] for key in df_container.keys()]) # first element of tuple is method
        self.XGs = set([key[1] for key in df_container.keys()]) # DG or SG
        self.additive_threshold = additive_threshold

    def fill_container(self):
        """Fill df_container with drug combination effect classification data frames for each method and XG in container dict."""
        for method in self.methods:
            dcd = DrugCombinationEffectClassification(
                self.df_container.get((method, 'SG'), None),
                self.df_container.get((method, 'DG'), None),
                self.additive_threshold) # each key level are dict of dfs
            dcd.get_gr_df()
            dcd.get_normalized_gr_df()
        return self.df_container
    
    def write_to_csv(self):
        """Save data frames in df_container in data directory."""
        def filename_format(data_directory, method, XG, df_type):
            """Format file name for saving data frames in data directory.
            File will also be used in data_analysis section for carbon allocation
            """
            if df_type == 'biomass':
                return os.path.join(data_directory, f'BM_{XG}_{method}.csv')
            if df_type == 'flux':
                return os.path.join(data_directory, f'flux_{XG}_{method}.csv')
            if df_type == 'growth_rate':
                return os.path.join(data_directory, f'gr_{XG}_{method}.csv')
            if df_type == 'normalized_growth_rate':
                return os.path.join(data_directory, f'normalized_gr_{XG}_{method}.csv')
            if df_type == 'drug_response_classification':
                return os.path.join(data_directory, f'drug_response_classification_{method}.csv')
            print('Unexpected df_type, check argument order')

        for df_type in ['biomass', 'growth_rate', 'flux', 'normalized_growth_rate', 'drug_response_classification']:
            for method, XG in itertools.product(self.methods, ['SG', 'DG']):
                df = self.df_container[method, XG].get(df_type, None) # SG do not have classification key
                if df is not None:
                    file_path = filename_format(self.data_directory, method, XG, df_type)
                    df.to_csv(file_path)
                    logger.info(f'Saved data frames for {method, XG}: biomass, flux, growth_rate, normalized_growth_rate, drug_response_classification in {file_path}')

        # write flux into one file
        for method in self.methods:
            flux_df_full = pd.concat([self.df_container[method, XG]['flux'] for XG in self.XGs])
            file_path = os.path.join(self.data_directory, f'flux_analysis_{method}.csv')
            flux_df_full.to_csv(file_path)
            print(f'biomass, growth_rate, normalized_growth_rate, drug_response_classification, flux data derived from {method} were saved in {file_path}')
