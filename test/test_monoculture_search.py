from scarcc.preparation.perturbation.alpha_finder import get_alpha_biomass_df
from scarcc.preparation.metabolic_model import BasicModel
import os
def find_directory(start_directory, directory_name):
    current_directory = os.path.abspath(start_directory)

    while True:
        directory_path = os.path.join(current_directory, directory_name)
        if os.path.isdir(directory_path):
            return directory_path

        # Move up one level
        current_directory = os.path.dirname(current_directory)

        # Check if reached the root directory
        if current_directory == os.path.dirname(current_directory):
            # Directory not found
            return None
data_directory = find_directory(os.path.abspath(__file__), 'models')
E0, S0, all_components = BasicModel(model_directory=data_directory, flux_weighting=True).load_ES_models()

df_list = get_alpha_biomass_df([S0], 0.5, potential_genes=['folA'])
print(df_list)




