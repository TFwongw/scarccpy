from setup import convert_arg_to_list
import itertools
from flux_snapshot import get_SG_DG

def plot_coculture(Biomass_df, gene_list,culture ='coculture', species_id_list = ['E0', 'S0']):
    species_id_list = convert_arg_to_list(species_id_list)

    if not isinstance(species_id_list[0], str):
        species_id_list = [m.id for m in species_id_list]
    colnames = list()
#     for model, current_gene in itertools.product([E0, S0],convert_arg_to_list(gene_list)):
    for model_id, current_gene in itertools.product(species_id_list,convert_arg_to_list(gene_list)):
        colnames.append('_'.join([model_id, current_gene, culture]))
    
    Biomass_df[colnames].plot()
    return None

# def plot_SG_DG_compare(Biomass_df, SG, DG_list, **kwargs):
#     SG_DG = get_SG_DG(DG_list, explode=True)
#     SG=convert_arg_to_list(SG)
#     DG = list(SG_DG.query('SG in @SG').index)
#     if DG:
#         print(SG, 'gcomb growth curve')
#         plot_coculture(Biomass_df, list(itertools.chain(*[DG,SG, ['Normal']])), **kwargs)
#     else: print(SG, 'do not have gcomb')

def plot_SG_DG_compare(Biomass_df, SG, DG_list, **kwargs):
    DG_list = convert_arg_to_list(DG_list)
    SG_DG = get_SG_DG(DG_list, explode=True)
    if SG is None:
        SG = list(itertools.chain(*[gcomb.split('.') for gcomb in DG_list]))
    else:
        SG=convert_arg_to_list(SG)

    DG = list(SG_DG.query('SG in @SG').index)
    if DG:
        print(list(itertools.chain(*[DG,SG, ['Normal']])))
        plot_coculture(Biomass_df, list(itertools.chain(*[DG,SG, ['Normal']])), **kwargs)
    else: print(SG, 'do not have gcomb')
# plot_coculture(Biomass_df, ['murA', 'murA.pyrD'],culture='monoculture')
# plot_SG_DG_compare(SG, culture='monoculture',species_list=['E0'])