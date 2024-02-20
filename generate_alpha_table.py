import os
# from scarcc.preparation.perturbation.alpha_finder import get_div_obj_df
# from scarcc.preparation.metabolic_model import BasicModel

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
            print(f'Directory {directory_name} not exist, make directory at root')
            os.makedirs(directory_name)

# data_directory = find_directory(os.path.abspath(__file__), 'models')
data_directory = find_directory(os.path.abspath(__file__), 'fodels')
# E0, S0, all_components = BasicModel(model_directory=data_directory, flux_weighting=True).load_ES_models()

# df_list = get_div_obj_df([S0], 0.4, potential_genes=['folA'])
# print(df_list)
# print('Script run without error')

# # make DATA directory
# file = os.path.abspath(__file__)
# data_directory = os.path.join(os.path.dirname(file), 'Data')


# os.makedirs('models', exist_ok=True)
print(os.path.abspath(__file__))
print(data_directory)
