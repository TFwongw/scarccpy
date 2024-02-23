"""Summary statistics from biomass growth record"""

import pandas as pd
from typing import List, Union
from scarcc.util import convert_arg_to_list
from ..flux.flux_snapshot import set_GI_SP_as_MI

def turn_to_normal_colname(colname):
    colname = colname.split('_')
    new_colname = '_'.join([colname[0], 'Normal', colname[2]])
    return new_colname

def convert_col_to_gene_inhibition(s):
    if len(s.split('.')) <=2:
        _, gcomb, _ = s.split('_')
        return gcomb
    _, gcomb, _, lv_pairs = s.split('_')
    return gcomb + '_' + lv_pairs

def get_biomass_df(file_path_list: Union[str, List]):
    """read biomass_df from file"""
    biomass_df =  pd.concat(
            [pd.read_csv(file, index_col='cycle')
            for file in convert_arg_to_list(file_path_list)]
        ,axis=1)
    if biomass_df.columns.str.contains('Normal').any():
        return biomass_df
    
    normal_columns = [ele for ele in biomass_df.columns if '0.0' in ele]
    normal_columns_rename_dict = {col: turn_to_normal_colname(col) for col in normal_columns}
    biomass_df.rename(columns=normal_columns_rename_dict, inplace=True)
    return biomass_df

def get_XG_cycle_from(desired_cycle):
    """separate SG and DG cycle from desired_cycle"""
    if len(desired_cycle.index[-1].split('.')) <=2:
        SG_cycle = desired_cycle.loc[[len(ele.split('.')) ==1 for ele in desired_cycle.index]]
        DG_cycle = desired_cycle.loc[[len(ele.split('.')) ==2 for ele in desired_cycle.index]]
    else:
        SG_cycle = desired_cycle.loc[['0' in ele for ele in desired_cycle.index]]
        DG_cycle = desired_cycle.loc[['0' not in ele for ele in desired_cycle.index]]
    SG_cycle.index.name='SG'
    DG_cycle.index.name='DG'
    SG_cycle.columns.name=None
    DG_cycle.columns.name=None
    return SG_cycle, DG_cycle

def search_gr_cycle_with_biomass(df_search, biomass_values):
    """search for cycle with biomass value
    biomass value at half maximum represents middle of exponential phase
    
    Parameters
    ----------
    df_search: pd.Series
        biomass record
    biomass_values: list
        biomass values to search

    Returns
    -------
    list of boolean if biomass record is larger than queried biomass value
    """
    return [df_search[df_search >= biomass_value].idxmin()
                for biomass_value in list(biomass_values)]

def get_maximum_growth_cycle(desired_biomass, scale_biomass_diff=0.1):
    """get summary statistics form biomass record
    
    Parameters
    ----------
    desired_biomass: pd.DataFrame
        biomass record
    scale_biomass_diff: float
        scale for biomass difference(max-min) to account for noise, 
        delayed start and advanced end cycle for exponential phase
        
    Returns
    -------
    c_max_gr: cycle at middle of exponential phase
    start: growth start cycle
    end: growth end cycle
    bool_growing: boolean if growing, cell not reach stationary phase
    """
    c_max_gr = desired_biomass.iloc[1]+ (desired_biomass.iloc[-1] - desired_biomass.iloc[1])/2
    bool_growing = ((desired_biomass.iloc[-1]-desired_biomass.iloc[-5])/desired_biomass.iloc[-1]).apply(lambda x: x > 1e-10)
    for k, bool_grow in bool_growing.items():
        if bool_grow:
            c_max_gr[k] = desired_biomass[k].iloc[-6]
    biomass_diff = (desired_biomass.iloc[-1]-desired_biomass.iloc[0])
    start = desired_biomass.iloc[0] + biomass_diff*scale_biomass_diff
    end = desired_biomass.iloc[0] + biomass_diff*(1-scale_biomass_diff)
    return c_max_gr, start, end, bool_growing

def get_desired_cycle(biomass_df, log_step=5, scale_biomass_diff=0.1):
    """get desired growth cycle from biomass record
    
    Parameters
    ----------
    biomass_df: pd.DataFrame
        biomass record
    log_step: int
        log step of flux record
    scale_biomass_diff: float
        argument for get_maximum_growth_cycle, refined exponential phase
        
    Returns
    -------
    desired_cycle: pd.DataFrame
        summary statistics from biomass record
    """
    def correct_cycle(cycle): # round to nearest available flux log time point
        if cycle < log_step:
            return log_step
        return round(cycle / log_step) * log_step

    def get_growth_phase_length():
        return ((desired_cycle['end'] - desired_cycle['start'])*(1-desired_cycle.bool_growing) + # if not growing, not changing growth phase length
                1e4*(desired_cycle.bool_growing)) #if growing, set growth length to 999
    
    def split_index_to_cols(df):
        items = df.index.str.split('_')
        columns = ['Species', 'Gene_inhibition', 'culture']
        if len(items[0]) > 3:
            columns.extend(['alpha_lv_pairs'])
        return pd.DataFrame(items.tolist(), index=df.index, columns=columns)
    desired_biomass_df = pd.DataFrame(get_maximum_growth_cycle(biomass_df, scale_biomass_diff=scale_biomass_diff)
                                    , index=['c_max_gr', 'start', 'end', 'bool_growing'])
    
    desired_cycle = (desired_biomass_df.iloc[:-1]
                .apply(lambda x: 
                        search_gr_cycle_with_biomass(biomass_df.loc[:,x.name],x))
                .T)
    desired_cycle['bool_growing'] = desired_biomass_df.T.bool_growing
    desired_cycle['cycle_max_gr'] = desired_cycle['c_max_gr'].apply(correct_cycle) # -> cycle_max_gr
    desired_cycle['growth_phase'] = desired_cycle[['start', 'end']].values.tolist()
    desired_cycle['growth_phase_length'] = get_growth_phase_length()
    desired_cycle['end_cycle'] = biomass_df.index[-1]//5*5
    desired_cycle = desired_cycle.join(split_index_to_cols(desired_cycle))
    if len(desired_cycle.Gene_inhibition.unique()) >1:
        desired_cycle = desired_cycle.set_index('Gene_inhibition')[['cycle_max_gr', 'bool_growing', 'growth_phase','growth_phase_length', 'Species','culture','end_cycle']]
    else:
        desired_cycle.index = ['_'.join([x[1],x[3]]) for x in desired_cycle.index.str.split('_')]
        desired_cycle.Gene_inhibition = desired_cycle.index

    return desired_cycle

def convert_col_to_gene_inhibition(biomass_column):
    if len(biomass_column.split('.')) <=2:
        _, gcomb, _ = biomass_column.split('_')
        return gcomb
    _, gcomb, _, lv_pairs = biomass_column.split('_')
    return gcomb + '_' + lv_pairs

def get_end_BM(biomass_df):    
    """retrieve biomass at end cycle"""
    def get_species_frac_binned(end_biomass): # derive species relative abundance in media
        stable_species_frac = end_biomass.reset_index('Species', drop=True)['BM_consortia_frac']
        # stable_species_frac = end_biomass.query('Species=="E0"').reset_index('Species', drop=True)['BM_consortia_frac'] # whyv need query E0

        labels = ['S', 'slight E', 'E']
        bins = [0, 0.5, 0.7, 1]
        stable_species_frac_binned = pd.cut(stable_species_frac, bins=bins, labels=labels)
        return stable_species_frac_binned
    
    def add_to_end_biomass(end_biomass): # additional columns counting species relative population 
        def get_ratio_and_std_col(model_id):
            temp_df = separate_Species_df(end_biomass, model_id, inc_Species=True)
            temp_df['BM_consortia_frac'] = temp_df.Biomass/temp_df.Total_BM
            temp_df['standardized_BM'] = temp_df.Biomass/temp_df.loc['Normal', 'Biomass']
            return temp_df
        
        def separate_Species_df(df, model_id, inc_Species = False):
            """separate species from df"""
            def get_species_loc(model_id):
                return list(temp_df.Species == model_id)
            
            if not isinstance(model_id, str):
                model_id = model_id.id
            df = df.loc[:,~df.columns.str.contains('Ndiff|ESdiff')]
            temp_df = pd.concat([df['Species'],df.select_dtypes(include = ['float'])], axis=1)
            result_df = temp_df.loc[get_species_loc(model_id)]
            return result_df if inc_Species else result_df.drop('Species',axis=1)
        return pd.concat([get_ratio_and_std_col('E0'), get_ratio_and_std_col('S0')]).set_index('Species', append=True).drop(['Total_BM', 'Biomass'], axis=1)

    biomass = pd.DataFrame(biomass_df.loc[:, [col for col in biomass_df.columns if 'coculture' in col]].iloc[-1])
    biomass.columns = ['Biomass']
    biomass['Species'] = list(pd.Series(biomass.index).apply(lambda x: str(x).split('_')[0]))
    biomass['Gene_inhibition'] = list(pd.Series(biomass.index).apply(convert_col_to_gene_inhibition))
    biomass = biomass.set_index(['Gene_inhibition'])
    biomass['Total_BM'] = biomass.groupby('Gene_inhibition').Biomass.sum()
    end_biomass = set_GI_SP_as_MI(biomass).join(set_GI_SP_as_MI(add_to_end_biomass(biomass)))
    binned_col = get_species_frac_binned(end_biomass)
    # return end_biomass, binned_col

    # end_biomass = end_biomass.merge(binned_col, left_index=True, right_index=True, suffixes=('', '_binned'))
    # end_biomass['BM_consortia_frac_binned'] = end_biomass['BM_consortia_frac_binned'].cat.add_categories(['No growth'])
    end_biomass['BM_consortia_frac_binned'] = list(binned_col)
    end_biomass.loc[end_biomass.Biomass<1.01e-8,'BM_consortia_frac_binned'] = 'No growth'
    return end_biomass
