import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from scarcc.utils import convert_arg_to_list

def scatter_xycol_A(flux_analysis, x_col='EX_bulk_ac_e', y_col='BIOMASS_iRR1083_metals'):
    sns.set_context('paper', font_scale=1.5)
    # cross_m = flux_analysis.query('XG=="SG"').query("culture=='coculture'")[['EX_bulk_ac_e', 'BIOMASS_iRR1083_metals']]
    cross_m = flux_analysis.set_index('Species', append=True)[['EX_bulk_ac_e', 'BIOMASS_iRR1083_metals', 'XG']]
    cross_E = cross_m.query('Species=="E0"')[[x_col]].reset_index('Species',drop=True)
    cross_S = cross_m.query('Species=="S0"')[[y_col]].reset_index('Species',drop=True)
    cf = cross_E.merge(cross_S, left_index=True, right_index=True).dropna()
    cf = cf.merge(flux_analysis[['XG']], left_index=True, right_index=True, how='left')

    sns.lmplot(data=cf, x=x_col, y=y_col, hue='XG')
    plt.title('Cross Feeding Metabolite Production Fluxes')
    plt.xlabel('E.coli Acetate Secretion')
    plt.ylabel('S. enterica Methionine Secretion')

def scatter_xycol_B(df, x_col, y_col, prefix='percent', XG=None):
    rename_label_dict = {'EX_bulk_ac_e': 'Acetate', 'EX_co2_e': 'CO2', 
                        'BIOMASS_Ec_iML1515_core_75p37M': 'Biomass', 'Waste': 'Waste'}
    color =dict(zip(['Antagonistic', 'Additive', 'Synergistic'], ['#25fe57', '#ebeb0f', '#fba2f7']))
    if XG is not None:
        XG=convert_arg_to_list(XG)
        df = df.query('XG in @XG')
    plt.figure(figsize=(9, 8))
    sns.scatterplot(data=df, x='_'.join([prefix, x_col]), y='_'.join([prefix, y_col]), 
                    hue='Drug_comb_effect_coc', palette=color, alpha=.8, s=100) 
    
    plt.title('Change in Normalized Carbon Flux \n (Double - Single Gene Inhibition)')
    # plt.xlabel('Change in Acetate production')
    # plt.ylabel('Change in CO2 production')
    plt.xlabel(f'Change in {rename_label_dict[x_col]} Production')
    plt.ylabel(f'Change in {rename_label_dict[y_col]} Production')
    plt.legend(prop={'size': 14}, loc='upper left', title='Drug Combination Effect', title_fontsize=14)
