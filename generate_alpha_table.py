# import pandas as pd
# from scarcc.preparation.metabolic_model import BasicModel
# from scarcc.preparation.find_directory import find_directory

from importlib import reload
import scarcc.preparation.perturbation.alpha_finder.monoculture
reload(scarcc.preparation.perturbation.alpha_finder.monoculture)

import scarcc.preparation.perturbation.alpha_finder.alpha_finder
reload(scarcc.preparation.perturbation.alpha_finder.alpha_finder)
from scarcc.preparation.perturbation.alpha_finder.alpha_finder import AlphaFinderConfig

from scarcc.preparation.perturbation.alpha_finder.monoculture import MonocultureAlphaFinder, get_div_obj_df
# data_directory = find_directory(os.path.abspath(__file__), 'models')

model_directory = find_directory('models', __file__)
# data_directory = find_directory('Data', __file__)
# alpha_table = pd.read_csv(os.path.join(data_directory, 'alpha_table_m1.csv'), index_col=0)
# E0, S0, all_components = BasicModel(model_directory=model_directory, flux_weighting=True).load_ES_models() 
df = get_div_obj_df([E0], target_obj_val=0.5, potential_genes=['folA'])
df
# AFC = AlphaFinderConfig(data_directory)

# model = 1
# search_alpha = 1
# current_gene = 'folA'
# target_obj_val = 0.5
# precision = 2
# acceptance_threshold_upper = 0.9
# acceptance_threshold_lower = 0.3
# # MAF = MonocultureAlphaFinder(model_list, 1, 'folA')

# AF = MonocultureAlphaFinder(model=model,
#             search_alpha = search_alpha,
#             current_gene = current_gene, 
#             target_obj_val = target_obj_val, 
#             exp_leap=3,
#             precision=precision,
#             acceptance_threshold_upper = acceptance_threshold_upper,
#             acceptance_threshold_lower = acceptance_threshold_lower)

