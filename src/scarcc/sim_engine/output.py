import pandas as pd
import os
import itertools
import logging

from scarcc.data_analysis import convert_po_col
from scarcc.data_analysis.growth.growth_rate import get_growth_rate_df

logger = logging.getLogger(__name__)

def merge_single_gene_gr(sg_df, dg_df):
    # gcomb as index, retrive the corresponding single gene in sg_df
    # ? Normal col && normal str correspondence in dict converion first for checker board?
    
    dg_df = pd.concat([sg_df.loc[['Normal']], dg_df]).rename(index={'Normal': 'Normal.Normal'})
    dg_df['First_gene'], dg_df['Second_gene'] = list(zip(*dg_df.index.str.split('.'))) # ? Need handel checkerboard, wirte into function
    dg_df.rename(index={'Normal.Normal': 'Normal'}, inplace=True)
    
    dg_df = (dg_df
                .merge(sg_df, left_on='First_gene', right_index=True, suffixes=['','_First_gene'])
                .merge(sg_df, left_on='Second_gene', right_index=True, suffixes=['','_Second_gene'])
                )
    return dg_df

def add_additive_and_drug_comb_response(df): # after normalization of DG_gr
    for sp_cul in itertools.product(['E0', 'S0'], ['coculture', 'monoculture']):
        sp_cul = '_'.join(sp_cul)
        df[f'Predicted_additive_effect_{sp_cul}'] = df[f'{sp_cul}_First_gene'] * df[f'{sp_cul}_Second_gene']
        df[f'po_diff_{sp_cul}'] = df[sp_cul] - df[f'Predicted_additive_effect_{sp_cul}'] # ! op_col
        df[f'Drug_comb_effect_{sp_cul}'] = convert_po_col(df[f'po_diff_{sp_cul}'], additive_threshold=0.05)
    return df

def reorder_columns(df):
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
    def __init__(self, dfs_in_SG_layer, dfs_in_DG_layer):
        self.dfs_in_SG_layer = dfs_in_SG_layer
        self.dfs_in_DG_layer = dfs_in_DG_layer
        self.generate_DG = dfs_in_DG_layer is not None
        self.pure_dg_gr = None

    # construct a sub her, normalize df
    def get_gr_df(self):
        sg_df = get_growth_rate_df(self.dfs_in_SG_layer['biomass'])
        self.dfs_in_SG_layer['growth_rate'] = reorder_columns(sg_df)
        if self.generate_DG:
            dg_df = get_growth_rate_df(self.dfs_in_DG_layer['biomass'])
            self.pure_dg_gr = dg_df.copy() # This is before merge, because joining sg_df in non-normalized and normalized 

            dg_df = merge_single_gene_gr(sg_df, dg_df)
            self.dfs_in_DG_layer['growth_rate'] = reorder_columns(dg_df)

    def get_normalized_gr_df(self):
        # sg_df = self.dfs_in_SG_layer['growth_rate'].copy() # reassignment fo not affect non-normalized gr
        sg_df = self.dfs_in_SG_layer['growth_rate']
        sg_df = sg_df.div(sg_df.loc['Normal'])
        self.dfs_in_SG_layer['normalized_growth_rate'] = reorder_columns(sg_df)
            
        if self.generate_DG:
            dg_df = self.pure_dg_gr
            dg_df = dg_df.div(sg_df.loc['Normal'])
            dg_df = merge_single_gene_gr(sg_df, dg_df)
            dg_df = add_additive_and_drug_comb_response(dg_df)
            self.dfs_in_DG_layer['normalized_growth_rate'] = reorder_columns(dg_df)
            self.dfs_in_DG_layer['drug_response_classification'] = dg_df.filter(regex='Drug_comb_effect_')

class MethodDataFiller:
    # TODO: handler for not supplying SG biomass, OR handle for the first construction of container dict
    # ! Is in fact op diff rather than po diff - rename column
    # for each method, cal the ~. , og idea go into the method, with SG and DG as sub dict
    # one level dict - query method as first element in tuple
    # ！if DG none then no cal any for DG OR ensure both SG and DG exist
    def __init__(self, df_container, data_directory): #? make a base accept only sg, dg
        self.df_container = df_container
        self.data_directory = data_directory
        self.methods = set([key[0] for key in df_container.keys()]) # first element of tuple is method
        self.XGs = set([key[1] for key in df_container.keys()]) # DG or SG

    def fill_container(self):
        for method in self.methods:
            dcd = DrugCombinationEffectClassification(self.df_container[method, 'SG'], self.df_container[method, 'DG']) # each key level are dict of dfs
            dcd.get_gr_df()
            dcd.get_normalized_gr_df()
        return self.df_container
    
    def write_to_csv(self):
        def filename_format(data_directory, method, XG, df_type):
            if df_type == 'biomass':
                return os.path.join(data_directory, f'BM_{XG}_{method}.csv')
            if df_type == 'growth_rate':
                return os.path.join(data_directory, f'gr_{XG}_{method}.csv')
            if df_type == 'normalized_growth_rate':
                return os.path.join(data_directory, f'normalized_gr_{XG}_{method}.csv')
            if df_type == 'drug_response_classification':
                return os.path.join(data_directory, f'drug_response_classification_{method}.csv')
            print('Unexpected df_type, check argument order')

        for df_type in ['biomass', 'growth_rate', 'normalized_growth_rate', 'drug_response_classification']:
            for method, XG in itertools.product(self.methods, ['SG', 'DG']):
                df = self.df_container[method, XG].get(df_type, None) # SG do not have classification key
                if df is not None:
                    file_path = filename_format(self.data_directory, method, XG, df_type)
                    df.to_csv(file_path)
                    logger.info(f'Saved data frames for {method, XG}: biomass, growth_rate, normalized_growth_rate, drug_response_classification in {file_path}')

        # write flux into one file
        for method in self.methods:
            flux_df_full = pd.concat([self.df_container[method, XG]['flux'] for XG in self.XGs])
            flux_df_full.to_csv(os.path.join(self.data_directory, f'flux_analysis_{method}.csv'))
            print(f'biomass, growth_rate, normalized_growth_rate, drug_response_classification, flux data derived from {method} were saved in {self.data_directory}')
