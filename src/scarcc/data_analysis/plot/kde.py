import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def get_fs_kde_plot(fs_change):
    drop_list = ['percent_EX_lcts_e', 'percent_EX_for_e', 'percent_EX_gal_e', 'percent_EX_glyclt_e', 'percent_EX_met__L_e','percent_EX_hacolipa_e']
    fs_plot = fs_change.drop([col for col in fs_change.columns if col in drop_list], axis=1)
    fs_plot = pd.melt(fs_plot, id_vars=['Drug_comb_effect_coc'], var_name='reaction', value_name='flux')
    return fs_plot

def plot_kde(plot_df, reaction=None, col_prefix='percent_', common_norm=False, ylim=None): 
    def display_label(k):
        if 'Waste' in k:
            return 'Waste'
        if 'co2' in k:
            return 'CO2'
        if 'ac_e' in k:
            return 'Acetate'
        if 'BIOMASS' in k:
            return 'Biomass'
    
    color_dict = {
    'Waste': 'purple',
    'CO2': 'green',
    'Acetate': 'blue',
    'Biomass': 'red',
    }
    # color_dict = {col_prefix+k: v for k, v in color_dict.items()}
    sns.set_context('paper', font_scale=1)
    plot_df = plot_df.query("reaction.str.contains(@col_prefix)").copy()
    plot_df.rename(columns={'Drug_comb_effect_coc': 'Drug Combination Effect'}, inplace=True)
    plot_df['reaction'] = plot_df['reaction'].map(display_label)
    normalized_flux = 'Normalized Carbon Flux'
    plot_df.rename({'flux': normalized_flux}, inplace=True, axis=1)
    print(plot_df.columns)
    g = sns.FacetGrid(plot_df, col='Drug Combination Effect')
    transparency=0.15
    g.map_dataframe(sns.kdeplot, x=normalized_flux, hue="reaction",
                    fill=True, common_norm=common_norm, alpha=transparency, palette=color_dict)
    # g.map_dataframe(sns.kdeplot, x="flux", hue="reaction", fill=True, common_norm=True, alpha=0.2, palette=color_dict)
    if ylim:
        g.set(ylim=ylim)
    title = reaction if reaction else 'Difference in Normalized Carbon Secretion Profile (Double Gene minus Single Gene)'
    g.fig.suptitle(title)
    
    patches = [mpatches.Patch(color=v, label=k) for k,v in color_dict.items()]
    leg = plt.legend(handles=patches)

    for lh in leg.legend_handles: 
        lh.set_alpha(0.25)
    plt.tight_layout()