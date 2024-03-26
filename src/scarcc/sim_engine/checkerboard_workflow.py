"""Checkerboard workflow for running simulation with checkerboard alpha table"""

import os
import logging
import pandas as pd
import concurrent.futures
import ast

from scarcc.sim_engine.output import MethodDataFiller
from .simulation_workflow import get_simulation_output, unpack_future_result_per_key, concat_result

logger = logging.getLogger(__name__)

def convert_checker_alpha_table_to_dict(alpha_table_checkerboard: pd.DataFrame): 
    """Convert checkerboard alpha table to dictionary using lv_pairs as key"""
    d= dict()
    if isinstance(alpha_table_checkerboard.lv_pairs.iloc[0], str):
        alpha_table_checkerboard['lv_pairs'] = alpha_table_checkerboard.lv_pairs.apply(ast.literal_eval)
    for lv_pairs in alpha_table_checkerboard.lv_pairs.unique():
        d['_'+'.'.join(map(str, lv_pairs))] = alpha_table_checkerboard.query('lv_pairs == @lv_pairs')
    return d

# TODO: import from SG data that already ran/ BUT IF result has only levels are recorded
def run_checkerboard_workflow(alpha_table: pd.DataFrame, data_directory: str,  max_cpus=12, additive_threshold=0.05, **kwargs):
    """Run simulation with checkerboard alpha table for 2 target genes
    
    Parameters
    ----------
    alpha_table : pd.DataFrame
        Alpha table with checkerboard pattern
    data_directory : str
        Directory to save the simulation results
    max_cpus : int
        Maximum number of cpus to use for simulation
    additive_threshold : float
        Threshold for drug combination response to be considered non-additive    
    """

    def added_Normal_cols(biomass_df):
        """Append Normal columns to the biomass data frame as control row"""
        def rename_to_Normal(colname):
            splitted = colname.split('_')
            splitted[1] = 'Normal'
            splitted[-1] = '0.0'
            return '_'.join(splitted)
        Normal_df = biomass_df[[col for col in biomass_df.columns if col.endswith('_0.0')]].copy()
        Normal_df.columns = Normal_df.columns.map(rename_to_Normal)
        return pd.concat([Normal_df, biomass_df], axis=1)

    method = 'checkerboard'
    available_cpus = min(os.cpu_count(), max_cpus)
    task_dict = {}
    current_gene = list(alpha_table.index.unique())
    if len(current_gene) > 2:
        raise ValueError('Table contains more than 2 genes, '
                        'please filter the appropriate gene rows' 
                        'OR check if Gene_inhibition is already set as index')

    sub_alpha_dict = convert_checker_alpha_table_to_dict(alpha_table)
    SG_lv_alpha_dict = {k: v for k, v in sub_alpha_dict.items() if '0' in k}
    DG_lv_alpha_dict = {k: v for k, v in sub_alpha_dict.items() if '0' not in k}
    with concurrent.futures.ProcessPoolExecutor(max_workers=available_cpus) as executor:
        for XG, lv_alpha_dict in zip(['SG', 'DG'], [SG_lv_alpha_dict, DG_lv_alpha_dict]):
            task_dict[method, XG] = [
                executor.submit(
                    get_simulation_output,
                    current_gene=current_gene,
                    alpha_table=single_lv_alpha_table,
                    checker_suffix=checker_suffix,
                    **kwargs
                )
                for checker_suffix, single_lv_alpha_table in lv_alpha_dict.items() # each level is with key as '_lv1.lv2' and value as alpha_table at that lv_pair
            ]

    df_container = {key: [future.result() for future in future_list] for key, future_list in task_dict.items()}
    df_container = {key: unpack_future_result_per_key(result_list) for key, result_list in df_container.items()} # list of [biomass, flux] into {biomass : biomass_list, flux : flux_list}
    df_container = {k: concat_result(sub_container) for k, sub_container in df_container.items()} # column-wise for biomass, row-wise for flux data frame concatenation
    df_container[method, 'SG']['biomass'] = added_Normal_cols(df_container[method, 'SG']['biomass'])
    
    mdf = MethodDataFiller(df_container, data_directory, additive_threshold)
    mdf.fill_container()
    mdf.write_to_csv()
    return df_container