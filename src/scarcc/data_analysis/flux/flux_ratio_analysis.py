import re
import pandas as pd
import os 

from scarcc.preparation.metabolic_model import get_component
import itertools
from .flux_snapshot import (get_XG_cycle_from, set_GI_SP_as_MI)

if 'E0' not in globals():
    E0, S0, all_components = None, None,None

def get_metabolite_summary(df, metabolite, index=0, top_n=10, model=E0, all_components=all_components, concat=False): # disregard factor, because of reversible
    # already separated into production_flux and consumption_flux 
    
    if isinstance(df, pd.DataFrame):
        if len(df)>1:
            solution = df.loc[index] if isinstance(index, str) else df.iloc[index]
        else: 
            solution = df.iloc[0]
    elif isinstance(df, pd.Series):
        solution = df
    if isinstance(metabolite, str):
        metabolite = get_component(model, metabolite, all_components)
    
    r_in_sol = [r for r in metabolite.reactions if r.id in solution.index]
    flux = pd.DataFrame( # code modified from cobra
        data=[
            (
                # r.id,
                solution[r.id],
                r.get_coefficient(metabolite.id),
            )
            for r in r_in_sol
        ],
        columns=["flux", "factor"],
        # columns=["reaction", "flux", "factor"],
        index=[r.id for r in r_in_sol],
    )
    # Scale fluxes by stoichiometric coefficient. (positive is production, negative is consumption)
    flux["flux"] *= flux["factor"]
    flux.index.name = "reaction"
    flux['Speices'] = model.id
    flux['Gene_inhibition'] = solution.name
    production_flux = flux.query("flux > 0").copy()
    consumption_flux = flux.query("flux < 0").copy()
    production, consumption = production_flux["flux"].dropna().abs(), consumption_flux["flux"].dropna().abs()
    production_flux["percent"] = (production / production.sum()).apply("{:.2%}".format)
    consumption_flux["percent"] = 1*(consumption / consumption.sum()).apply("-{:.2%}".format) # negative sign as consumption

    # return consumption_flux.sort_values('flux', ascending=True)
    out_list = (production_flux.sort_values('flux', ascending=False)[:top_n],
            consumption_flux.sort_values('flux', ascending=True)[:top_n])
    return pd.concat(out_list) if concat else out_list

def get_antagonistic_df(syn_df):
    antagonistic_list = syn_df.loc[syn_df['P_O']<0].query("P_O<-0.01").index 
    # ? func get pwy col
    potential_pwy = (gene_combo_pathway
                    .loc[antagonistic_list,'Pathway']
                    .str.split(' \+ ') # series string need \+
                    .apply(lambda x: sorted(x))) 
    return pd.DataFrame(potential_pwy)

def remove_nan_from(x: pd.Series):
    return x.apply(lambda x: sorted(list(itertools.compress(x,[ele not in [None, np.nan] for ele in x]))))

def get_syn_df(gr_XG, model_id, gr_DG, additive_threshold = .01):
    diff_bins = [-10,-1*additive_threshold,1*additive_threshold,10]
    gr_bins = [-1,.2,1,10]
    
    if gr_XG.index.name == 'gene_inhibition':
        gr_XG.index.name = 'Gene_inhibition'
    syn_df  = pd.DataFrame([], index=gr_DG.index)
    syn_df['Predicted_growth_rate'] = gr_XG[f'Predicted_additive_effect_{model_id}_coculture']
    syn_df['Observed_growth_rate'] = gr_XG[f'{model_id}_coculture']
    syn_df['P_O'] = syn_df.Predicted_growth_rate - syn_df.Observed_growth_rate
    syn_df['Drug_comb_effect'] = pd.cut(syn_df['P_O'], bins=diff_bins, labels=['Antagonistic', 'Additive', 'Synergistic'])
    
    syn_df['PGR_bin'] = pd.cut(syn_df['Predicted_growth_rate'], bins=gr_bins, labels=['Low', 'Normal', 'High'])
    syn_df['OGR_bin'] = pd.cut(syn_df['Observed_growth_rate'], bins=gr_bins, labels=['Low', 'Normal', 'High'])
    
    syn_df.loc[(syn_df.Predicted_growth_rate < 1.5e-8) & (syn_df.Observed_growth_rate < 1.5e-8),'P_O'] = 0
    syn_df['gene_sort'] = (syn_df['P_O'] < 0)*10 + abs(syn_df['P_O'])
    syn_df['Species'] = model_id
    return syn_df

def get_p_o_df(gr_DG):
    p_o = pd.concat([get_syn_df(gr_DG, 'E0', gr_DG), get_syn_df(gr_DG, 'S0', gr_DG)])
    return p_o

# functions for complete flux df
def add_SG_to_p_o(desired_cycle, p_o_full):
    SG_cycle, _ = get_XG_cycle_from(desired_cycle)
    SG_list = list(SG_cycle.index.unique())
    if 'Normal' in SG_list:
        SG_list.remove('Normal')
    SG_empty = (pd.DataFrame(index=SG_list, columns=p_o_full.columns))
    SG_empty['Species'] = [['E0', 'S0']]*len(SG_empty)
    SG_empty = SG_empty.explode(['Species'])
    return pd.concat([p_o_full, SG_empty])

def get_full_df(desired_cycle, p_o_full, additional_df_list, Species='E0'): # Missing metabolite
#     df_list = [add_SG_to_p_o(), desired_cycle, end_BM, pwy_rct_df, flux_compare_df]
    p_o_w_SG = add_SG_to_p_o(desired_cycle, p_o_full)
    # additional_df_list = end_BM, pwy_rxn_df, flux_compare_df
    df_list = [p_o_w_SG, desired_cycle.query('culture=="coculture"')]
    df_list.extend(additional_df_list)
    merged_df = df_list.pop(0)
    for i, next_df in enumerate(df_list):
        manual_SI = False

        if 'Species' not in next_df:
            manual_SI = True
        if 'culture' in next_df:
            print(set(next_df.culture))
            next_df = next_df.query('culture=="coculture"')
            print(set(next_df.culture))
            if 'culture' in merged_df:
                next_df = next_df.drop('culture', axis=1) 
        if 'Species' in next_df:
            next_df = next_df.query('Species== @Species')
        merged_df = (set_GI_SP_as_MI(merged_df)
                    .join(set_GI_SP_as_MI(next_df), how='left')) # left join->DG only
        # inner_join SG info from desired cycle & alpha & flux
    return merged_df.reset_index(level='Species')