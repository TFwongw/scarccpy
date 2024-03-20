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
    def attach_lv_to_gene_inhibition(col_elements):
        Species, Gene_inhibition, Culture, lv_pair = col_elements
        return (Species, '_'.join([Gene_inhibition, lv_pair]), Culture)

    # derive growth rate for each column then retrieve the gene_inhibition as index
    # Example: Column E0_folA_coculture in Biomass_df column generated the growth rate and stored in column E0_coculture with index folA 
    df = Biomass_df.apply(get_growth_rate, axis=0).to_frame()
    df.columns = ['growth_rate']
    
    new_index = [i.split('_') for i in df.index]
    if len(new_index[-1]) == 4:
        new_index = [attach_lv_to_gene_inhibition(eles) for eles in new_index]
    index_names=['Species', 'Gene_inhibition', 'Culture']
    df.index = pd.MultiIndex.from_tuples(new_index, names=index_names)

    df = df.reset_index(['Species', 'Culture']).pivot(columns=['Species', 'Culture'], values='growth_rate')
    df.columns = ['_'.join(col) for col in df.columns]
    return df