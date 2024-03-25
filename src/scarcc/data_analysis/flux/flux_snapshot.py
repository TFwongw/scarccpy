import pandas as pd
from scarcc.preparation.metabolic_model import get_gene_id

def get_rcts_list(model, gcomb_list):
    """Get unique reactions list for each gene combination"""
    rcts_list = list()
    rcts_set = set()
    for i, gene in enumerate(gcomb_list):
        gene_rcts = [rct.id for rct in model.genes.get_by_id(get_gene_id(model, gene)).reactions]
        if i > 0:
            gene_rcts = list(set(gene_rcts) - rcts_set)
        rcts_list.append(gene_rcts)
        rcts_set = rcts_set | set(gene_rcts)
    return rcts_list

def adjust_flux_df(model, df, gene_combo: list, alpha_table:pd.DataFrame): # model for query of reactions
    """Adjust flux df with alpha value for each target gene in the gene combination.
    
    Adjustment for larger alpha come first and scale only once
    """
    # Don't use v1 cols to indicate gene_inhibition, will skip unidirectional reaction
    def query_alpha(gene_combo, alpha_table):
        splitted = gene_combo.split('.')
        is_checkerboard = len(splitted) > 2
        if not is_checkerboard:
            gcomb_alpha = {gene: alpha_table.loc[gene, f'{model.id}'] for gene in gene_combo.split('.')}
            return gcomb_alpha, alpha_table

        else:
            splitted = gene_combo.split('_')
            Gene_inhibition, lv_pair = [ele.split('.') for ele in gene_combo.split('_')]
            lv_pair = tuple([int(ele) for ele in lv_pair])
            gcomb_alpha = dict()
            alpha_table = alpha_table.query('lv_pairs == @lv_pair')
            for current_gene in Gene_inhibition:
                gcomb_alpha.update({current_gene:
                                    alpha_table.loc[current_gene, f'{model.id}']})
            return gcomb_alpha, alpha_table

    if 'Normal' not in gene_combo:
        v1_cols = df.filter(regex='v1').columns
        gcomb_alpha, alpha_table = query_alpha(gene_combo, alpha_table)
        gcomb_alpha = dict(sorted(gcomb_alpha.items(), key=lambda item: item[1], reverse=True))
        rcts_list = get_rcts_list(model, gcomb_alpha.keys()) # exclude repeated rct for gene with lower alpha

        for gene, rcts in zip(gcomb_alpha.keys(), rcts_list):
            for orig_col in rcts:
                alpha = alpha_table.loc[f'{gene}', f'{model.id}']
                v1_col = orig_col + "_v1"
                reversible = v1_col  in v1_cols
                if reversible:
                    df[f'{orig_col}'] = (df[f'{orig_col}'] + df[f'{v1_col}'])/alpha # only forward or backward != 0  
                    df = df.drop(f'{v1_col}', axis=1)  
                else:
                    df[f'{orig_col}'] = (df[f'{orig_col}'])/alpha # only forward or backward != 0              
    return df

def get_XG_cycle_from(desired_cycle):
    """Separate SG and DG for desired_cycle df."""
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

def set_GI_SP_as_MI(end_BM, manual_SI = False, multi_index = None):
    """Set Gene_inhibition and Species as MultiIndex, if they are in column or as index.
    
    If Species is not in columns, set Gene_inhibition as index.
    """
    if type(end_BM.index[0]) is not int: # Gene_inhibition is not index # Series.dtype -> dtype('int64')
        if end_BM.index.name in [None, 'DG']:
            end_BM.index.name = 'Gene_inhibition'
        if end_BM.index.name == 'Gene_inhibition':
            if 'Gene_inhibition' in end_BM.columns:
                end_BM = end_BM.drop('Gene_inhibition', axis=1)
            end_BM = end_BM.reset_index()
    if multi_index:
        return end_BM.set_index(multi_index)
    if 'Species' in end_BM.columns and manual_SI == False:
        return end_BM.set_index(['Gene_inhibition', 'Species'])
    else:
        return end_BM.set_index('Gene_inhibition')

def join_dfs_using_MI(df_list, how='left', multi_index = None):
    """Join dfs using MultiIndex."""
    df_list = [set_GI_SP_as_MI(df, multi_index=multi_index) for df in df_list]
    result_df = df_list[0]
    for df in df_list[1:]:
        result_df = result_df.join(df, how=how)
    return result_df

def get_SG_DG(DG_list, explode=True) -> pd.DataFrame:
    """create data frame with SG and DG columns from DG_list."""
    df = pd.DataFrame(pd.DataFrame(2*[DG_list], index=['Gene_inhibition','SG']).T.set_index('Gene_inhibition').SG.str.split('.'))
    return df.explode('SG') if explode else df
