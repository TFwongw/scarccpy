import pandas as pd
import re
from scarcc.preparation.metabolic_model import get_component
from scarcc.data_analysis import convert_po_col
# from scarcc.data_analysis.

def get_total_carbon(rxn, model, all_components, is_substrate=True,  carbon_detail=False, atp_only=False, atp_exclude=True): # TODO: exclude EX_met__L_e
    def get_n_carbon(formula):
        if formula in ['CO2', 'CH1O2']:
            return 1
        value = re.findall(r'C(\d+)', formula)
        numeric_value = int(value[0]) if value else 0
        return(numeric_value)
    
    def eval_include_in_dict(metabolite, factor, is_substrate=True):
        if atp_exclude and metabolite.id == 'atp_c':
            return False
        if rxn.id == 'BIOMASS_Ec_iML1515_core_75p37M' and atp_only and metabolite.id != 'atp_c':
            return False
        return (factor<0) == is_substrate and get_n_carbon(metabolite.formula)!=0
    
    if isinstance(rxn, str):
        rxn = get_component(model, rxn, all_components)
    d = {m.id: {'formula': m.formula,
                'formula_carbon': get_n_carbon(m.formula),
                'factor': factor,
                'carbon': factor * get_n_carbon(m.formula),}
                for m, factor in rxn.metabolites.items() if eval_include_in_dict(m, factor, is_substrate)}
    carbon_only = {k:v['carbon'] for k,v in d.items()}
    if carbon_detail:
        return carbon_only
    carbon_total = sum(carbon_only.values())
    return carbon_total

def get_carbon_allocation_summary(s, model, all_components, detailed=False, carbon_dict= {}, format=True): # accept series as input and return allocation df   
    for rxn, flux_value in s.filter(regex='EX_.*e$|BIOMASS_Ec_iML1515_core_75p37M$').dropna().items():
        
        carbon_per_flux = get_total_carbon(rxn, is_substrate=True, model=model, all_components=all_components, atp_exclude=True)

        carbon_dict.update({rxn: {
                            'carbon_per_flux': carbon_per_flux,
                            'flux_quantity': flux_value,
                            'carbon_exchange': float(flux_value * carbon_per_flux),
                            }})
    if detailed:
        return pd.DataFrame(carbon_dict).T

    carbon_allocation = pd.DataFrame({rxn: {'total_carbon': v['carbon_exchange']}
                        for rxn, v in carbon_dict.items() if abs(v['carbon_exchange'])>1e-4}).T
    normalize_carbon = carbon_allocation.loc['EX_lcts_e','total_carbon']
    carbon_allocation['percent'] = (carbon_allocation.total_carbon/normalize_carbon)
    waste_portion = -1-carbon_allocation.query('percent<0').percent.sum()
    if abs(waste_portion)<0.01: # if waste portion is too small, ignore it
        waste_portion = 0

    waste_row = pd.DataFrame([-1*waste_portion*normalize_carbon, waste_portion] # inclusion of waste product
                            ,columns=['Waste'], index=carbon_allocation.columns).T
    
    carbon_allocation = pd.concat([carbon_allocation, waste_row])
    if format:
        carbon_allocation['percent'] = carbon_allocation['percent'].apply(lambda x: f'{x:.2%}')
    carbon_allocation['Gene_inhibition'] = s.name if isinstance(s.name, str) else s.name[0]
    carbon_allocation.index.name = 'reaction'
    return carbon_allocation

def get_carbon_allocation_E_wide(E0, all_components, additive_threshold=0.05, flux_analysis_file='./Data/flux_analysis_m1.csv', gr_file='./Data/gr_DG_m1_normalized.csv', end_BM=None):
    # flux_analysis_full = pd.read_csv('./Data/flux_analysis_m1.csv', index_col=0)
    flux_analysis_full = pd.read_csv(flux_analysis_file, index_col=0)
    co_E = flux_analysis_full.query('Species=="E0"')
    if 'culture' in co_E.columns:
        co_E = co_E.query('culture=="coculture"')
    # if 'XG'in co_E.columns:
    #     co_E = co_E.query('XG=="DG"')
    gr_df = pd.read_csv(gr_file, index_col=0)
    # gr_df = pd.read_csv('./Data/gr_DG_m1_normalized.csv', index_col=0)
    carbon_allocation_E = pd.concat(co_E.apply(lambda x: 
                                            get_carbon_allocation_summary(x, E0, all_components,format=False), axis=1)
                                                .tolist())
    carbon_allocation_rxns = get_carbon_allocation_summary(co_E.iloc[0], E0, all_components).index # carbon allocation fluxes

    carbon_allocation_E_wide = carbon_allocation_E.reset_index().pivot(index='Gene_inhibition',columns='reaction',values=['total_carbon', 'percent'])
    # carbon_allocation_E_pivot = carbon_allocation_E.pivot(index='Gene_inhibition',columns='reaction',values=['total_carbon', 'percent'])
    carbon_allocation_E_wide.columns = ['_'.join([col[0], col[1]]) for col in carbon_allocation_E_wide.columns]
    carbon_allocation_E_wide = carbon_allocation_E_wide.merge(co_E[['XG']], left_index=True, right_index=True) # only E0
    partial_convert_po_col = lambda x: convert_po_col(x, additive_threshold=additive_threshold)
    carbon_allocation_E_wide['Drug_comb_effect_coc'] = partial_convert_po_col(gr_df['po_diff_E0_coculture'])
    carbon_allocation_E_wide['Drug_comb_effect_Emono'] = partial_convert_po_col(gr_df['po_diff_E0_monoculture'])
    carbon_allocation_E_wide['Drug_comb_effect_Smono'] = partial_convert_po_col(gr_df['po_diff_S0_monoculture'])

    # return carbon_allocation_E_wide
    if 'BM_consortia_frac_binned' not in carbon_allocation_E_wide.columns and end_BM is not None:
        carbon_allocation_E_wide = carbon_allocation_E_wide.merge(end_BM.query('Species=="E0"').BM_consortia_frac_binned, left_index=True, right_index=True)
    carbon_allocation_E_wide = carbon_allocation_E_wide.reset_index('Species', drop=True)
    # carbon_allocation_E_wide.to_csv('./Data/carbon_allocation.csv')
    return carbon_allocation_E_wide

## Difference in normalized carbon flux for double gene compared to First and second 
def get_percent_cols(df, col_suffix=None, col_prefix='percent'):
    df = df.filter(regex=col_prefix)
    if not col_suffix:
        return df
    df = df.filter(like=col_suffix) 
    df.columns = df.columns.str.replace(col_suffix,'')
    return df

def get_fs_change(gr_path, carbon_allocation_E_wide, additive_threshold=0.05):
    def get_sub_fs_change(fs_df, col_suffix='_First', col_prefix='percent'):
        fs = get_percent_cols(fs_df, col_suffix, col_prefix=col_prefix)
        # return carbon_allocation_E_wide
        fs_change = fs - get_percent_cols(carbon_allocation_E_wide.loc[fs.index], col_suffix=None, col_prefix=col_prefix) # change compared to single gene
        # fs_change = fs_change.merge(carbon_allocation_E_wide[['Drug_comb_effect_coc']], left_index=True, right_index=True).dropna(axis=0, how='all')
        return fs_change
    
    gr_df = pd.read_csv(gr_path, index_col=0)
    fs_df = carbon_allocation_E_wide.merge(gr_df[['First_gene', 'Second_gene']], left_index=True, right_index=True)

    fs_df = (fs_df.merge(carbon_allocation_E_wide, left_on='First_gene', right_index=True, suffixes=['','_First']) # add corresponding Frist and Second flux data
            .merge(carbon_allocation_E_wide, left_on='Second_gene', right_index=True, suffixes=['','_Second']))

    fs_change = pd.concat([get_sub_fs_change(fs_df, '_First'), get_sub_fs_change(fs_df, '_Second')])
    fs_change = fs_change.sort_index()
    fs_change['Nth_gene'] = list(['First','Second']*len(fs_change.index.unique()))
    # if 'Drug_comb_effect_coc' in fs_change.columns:
    #     fs_change = fs_change.drop(['Drug_comb_effect_coc'], axis=1)
    fs_change = fs_change.merge(carbon_allocation_E_wide.loc[:,'Drug_comb_effect_coc':], left_index=True, right_index=True)
    
    return fs_change