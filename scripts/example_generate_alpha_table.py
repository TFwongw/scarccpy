import os

from scarcc.preparation.metabolic_model import BasicModel
from scarcc.preparation.find_directory import find_directory
from scarcc.preparation.target_gene.gene_format_handler import get_DG_list, get_SG_list
from scarcc.preparation.perturbation.alpha_finder.monoculture import get_alpha_biomass_df

def main():
    # get file directory
    model_directory = find_directory('models', os.path.abspath(''))
    data_directory = find_directory('Data', os.path.abspath(''))

    # initialize model
    E0, S0, all_components = BasicModel(model_directory=model_directory, flux_weighting=True).load_ES_models()
    DG_list = get_DG_list(os.path.join(data_directory, 'GeneCombos.csv'), n_combos=None)
    SG_list = get_SG_list(DG_list)

    maf_kwargs = {
                'target_normalized_biomass': 0.5,
                'potential_genes': SG_list,
                'precision': 2,
                'acceptance_threshold_lower': 1}
    # default detailed_alpha_table=False
    alpha_table = get_alpha_biomass_df(model_list = [E0, S0], data_directory=data_directory, **maf_kwargs)
    print(alpha_table.head())

if __name__ == '__main__':
    main()