"""Utility functions for data processing and analysis."""

import pandas as pd
import re

def convert_arg_to_list(arg):
    """Converts argument to list if it is not already a list, tuple or set.
    
    Parameters
    ----------
    arg: str, list, tuple, set, pd.Series, pd.Index
    
    Returns
    -------
    list[arg] or arg in original format
    """
    if isinstance(arg, (pd.Series, pd.Index)):
        arg = list(arg)
    elif type(arg) not in [list, tuple, set]:
        arg = [arg] # differ from list(arg) -> conversion of str
    return arg

def rename_columns(df):
    """Rename columns of dataframe to remove special characters and spaces.
    Original column names are in the format: {Species}_{culture_characteristics}_{'[gene1, gene2]'}_{culture}
    renames to: {Species}_{gene1.gene2}_{culture}

    Parameters
    ----------
    df: pandas.DataFrame
    
    Returns
    -------
    pandas.DataFrame.columns    
    """
    df.columns = [re.sub('S0_ac_','S0.ac_', ele) for ele in df] # S0_ac -> S0.ac
    df.columns = [re.sub('S0_gal_','S0.gal_', ele) for ele in df] # S0_gal -> S0.ac
    df.columns = [re.sub(',','.',
           re.sub('\'|\(|\)| |\[|\]','',ele)) # ('gene1', 'gene2') -> gene1.gene2
           for ele in df.columns]
    return(df.columns)

def remove_Zero_col(df): # extend N differ than 0
    """Remove columns with all zero entries"""
    return(df.loc[:, ((df !=0) & (df.notnull())).any(axis=0)]) # ignore NA entry
