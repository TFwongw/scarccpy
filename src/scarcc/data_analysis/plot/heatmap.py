import pandas as pd
import numpy as np

def customize_row_color(df, col, color, col_elements=None):
    if 'comb' in col:
        lut = dict(zip(['Antagonistic', 'Additive', 'Synergistic'], ['#90fda9', '#ffff84', '#fed0fc']))
    else:
        elements = df[col].unique() if not col_elements else col_elements
        lut = dict(zip(elements, color))
    row_colors = df[col].map(lut)
    return row_colors

def generate_row_colors(fs_sub_change, color=[ '#4622fc','#9d22fc', '#da52d5', 'grey'], colnames=None, strip_DCE=True):
    row_colors_response = customize_row_color(fs_sub_change, 'Drug_comb_effect_coc', 'rbg')
    row_colors_Smono = customize_row_color(fs_sub_change, 'Drug_comb_effect_Smono', 'rbg')
    row_colors_Emono = customize_row_color(fs_sub_change, 'Drug_comb_effect_Emono', 'rbg')
    # row_colors_BM_frac = customize_row_color(fs_sub_change, 'BM_consortia_frac_binned',col_elements=['E', 'slight E','S'], color=['#ff05e8', '#e03fd8', '#f075e6'])
    # row_colors_BM_frac = customize_row_color(fs_sub_change, 'BM_consortia_frac_binned',col_elements=['E', 'slight E','S'], color=['#c14a09', '#D46C4E', '#F9AD6A'])
    row_colors_BM_frac = customize_row_color(fs_sub_change, 'BM_consortia_frac_binned',col_elements=['E', 'slight E','S','No growth'], color=color) # TODO: verify no growth adjustment 
    row_colors = pd.concat([row_colors_Smono,row_colors_response,row_colors_Emono, row_colors_BM_frac],axis=1)
    if colnames is None:
        row_colors.columns = ['DCE mono S', 'DCE co E', 'DCE monoc E', 'Rel. Abund.']
        if strip_DCE:
            row_colors.columns = ['mono S', 'co E', 'monoc E', 'Rel. Abund.']
    return row_colors

def relabel_xtick(clustered):
    xtick_relabel = {'Waste' : 'Waste',
                    'BIOMASS_Ec_iML1515_core_75p37M': 'Biomass',
                    'EX_bulk_ac_e': 'Acetate',
                    'EX_co2_e': 'CO2'}
    xtick_to_rename = [tx.get_text().replace('percent_','').replace('total_carbon_', '').replace('plot_', '') for tx in clustered.ax_heatmap.get_xticklabels()]
    xtick_to_rename = [xtick_relabel.get(ele, ele) for ele in xtick_to_rename]
    return xtick_to_rename

def relabel_clustermap(clustered, x_label_middle_str=''):
    clustered.ax_heatmap.set_yticklabels(clustered.ax_heatmap.get_yticklabels(), fontsize=8)
    clustered.ax_heatmap.set_ylabel('Double - Single Drug Inhibition', weight='bold', fontsize=12)
    clustered.ax_heatmap.set_xticklabels(relabel_xtick(clustered), rotation=0) # need check original label

    if x_label_middle_str:
        x_label_middle_str= x_label_middle_str+'\n'
    clustered.ax_heatmap.set_xlabel(f'''Difference in Normalized Carbon Flux \n{x_label_middle_str}(E. coli, Double - Single Drug Inhibition)''', weight='bold')
    return None

def assign_plot_total_E_wide(carbon_allocation_E_wide, log_waste=True):
    carbon_allocation_E_wide['plot_total_carbon_Waste'] = np.log(carbon_allocation_E_wide['total_carbon_Waste']) if log_waste else carbon_allocation_E_wide['total_carbon_Waste']
    carbon_allocation_E_wide.loc['Normal','plot_total_carbon_Waste']=-10

    # carbon_allocation_E_wide['plot_total_carbon_EX_bulk_ac_e'] = carbon_allocation_E_wide['total_carbon_EX_bulk_ac_e']/-2
    carbon_allocation_E_wide['plot_total_carbon_EX_bulk_ac_e'] = carbon_allocation_E_wide['total_carbon_EX_bulk_ac_e']/-1
    carbon_allocation_E_wide['plot_total_carbon_BIOMASS_Ec_iML1515_core_75p37M'] = carbon_allocation_E_wide['total_carbon_BIOMASS_Ec_iML1515_core_75p37M']*-1
    carbon_allocation_E_wide['plot_total_carbon_EX_co2_e'] = carbon_allocation_E_wide['total_carbon_EX_co2_e']*-1
    response_cols = ['Drug_comb_effect_coc', 'Drug_comb_effect_Emono', 'Drug_comb_effect_Smono', 'BM_consortia_frac_binned']
    keep_cols = ['plot_total_carbon_Waste', 'plot_total_carbon_EX_bulk_ac_e', 'plot_total_carbon_BIOMASS_Ec_iML1515_core_75p37M', 'plot_total_carbon_EX_co2_e']
    keep_cols.extend(response_cols)
    carbon_allocation_E_wide = carbon_allocation_E_wide[keep_cols]
    return carbon_allocation_E_wide