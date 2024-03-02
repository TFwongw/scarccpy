# accepts biomass as data frame then returns growth rate, normalized growth rate, classification of gene  combination effect
import os
import numpy as np
import pandas as pd
import itertools

from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

from scarcc.data_analysis import convert_po_col

def loglinear(x, a, b):
    return a*np.exp(b*(x))

def truncate_max_biomass(biomas_s):
    fit_col = biomas_s.dropna()
    cycles_to_fit = fit_col<fit_col.iloc[-1]*.9
    if any(cycles_to_fit.values) == False:
        return None
    
    truncate_col = fit_col.loc[cycles_to_fit]
    return truncate_col

def get_growth_rate(biomass_s, fit=False, plot_fitted=False, plot_og=False):
    """get growth rate from biomass series"""
    def get_fitted_curve(x, *popt):
        y_fit = loglinear(x, *popt)
        return y_fit
    
    def plot_curves(plot_fitted, plot_og=True, **kwargs):
        popt, pcov = curve_fit(loglinear, x, y, p0=[est_initial_pop, np_gr])

        if plot_og:
            plt.plot(x, y, label='Original Curve', **kwargs)
        if plot_fitted:        
            y_fit = get_fitted_curve(x,*popt)
            plt.plot(x, y_fit, label='Fitted Curve', **kwargs)
        
    fit_col = biomass_s.dropna()
    fit_col = truncate_max_biomass(biomass_s)

    if fit_col is None:
        return 0
    x = fit_col.index
    y= fit_col
    np_gr, est_initial_pop = np.polyfit(x, np.log(y), 1) # [B,A] -- log(y) = A+Bx
    plot_curves(plot_fitted, plot_og)
    return np_gr

def get_growth_rate_df(Biomass_df):
    # derive growth rate for each column then retrieve the gene_inhibition as index
    # Example: Column E0_folA_coculture in Biomass_df column generated the growth rate and stored in column E0_coculture with index folA 
    df = Biomass_df.apply(get_growth_rate, axis=0).to_frame()
    # df = pd.DataFrame(df).reset_index()
    df.columns = ['growth_rate']
    
    df.index = pd.MultiIndex.from_tuples([i.split('_') for i in df.index], names=['Species', 'Gene_inhibition', 'Culture'])
    df = df.reset_index(['Species', 'Culture']).pivot(columns=['Species', 'Culture'], values='growth_rate')
    df.columns = ['_'.join(col) for col in df.columns]
    return df

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
        df[f'po_diff_{sp_cul}'] = df[sp_cul] - df[f'Predicted_additive_effect_{sp_cul}']
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

class MethodDataFiller:
    # TODO: handler for not supplying SG biomass, OR handle for the first construction of container dict
    def __init__(self, df_container, data_directory):
        self.df_container = df_container
        self.data_directory = data_directory
        self.method = list(df_container.keys())[0]
        self.dfs_in_SG_layer = df_container[self.method]['SG']
        self.dfs_in_DG_layer = df_container[self.method]['DG']
        
    def get_gr_df(self):
        sg_df = get_growth_rate_df(self.dfs_in_SG_layer['biomass'])
        dg_df = get_growth_rate_df(self.dfs_in_DG_layer['biomass'])
        self.pure_dg_gr = dg_df.copy()

        dg_df = merge_single_gene_gr(sg_df, dg_df)        
        self.dfs_in_SG_layer['growth_rate'], self.dfs_in_DG_layer['growth_rate'] = [reorder_columns(df) for df in [sg_df, dg_df]]

    def get_normalized_gr_df(self):
        sg_df = self.dfs_in_SG_layer['growth_rate'].copy()
        dg_df = self.pure_dg_gr
        normal_row = sg_df.loc['Normal']
        sg_df, dg_df = sg_df.div(normal_row), dg_df.div(normal_row)
        dg_df = merge_single_gene_gr(sg_df, dg_df)
        dg_df = add_additive_and_drug_comb_response(dg_df)
        self.dfs_in_SG_layer['normalized_growth_rate'], self.dfs_in_DG_layer['normalized_growth_rate'] = [reorder_columns(df) for df in [sg_df, dg_df]]
        self.dfs_in_DG_layer['drug_response_classification'] = dg_df.filter(regex='Drug_comb_effect_')

    def fill_container(self):
        self.get_gr_df()
        self.get_normalized_gr_df()
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
            for XG in ['SG', 'DG']:
                df = self.df_container[self.method][XG].get(df_type, None)
                if df is not None:
                    file_path = filename_format(self.data_directory, self.method, XG, df_type)
                    df.to_csv(file_path)

        # write flux into one file
        flux_df_full = pd.concat([self.df_container[self.method][XG]['flux'] for XG in ['SG', 'DG']])
        flux_df_full.to_csv(os.path.join(self.data_directory, f'flux_analysis_{self.method}.csv'))
        print(f'All data frames derived from {self.method} were saved in {self.data_directory}')
