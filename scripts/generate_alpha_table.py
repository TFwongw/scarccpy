import os
import pandas as pd

from scarcc.preparation.metabolic_model import BasicModel
from scarcc.preparation.find_directory import find_directory 
from scarcc.preparation.perturbation.alpha_finder.monoculture import get_div_obj_df

# read potential genes
potential_genes = ['folA']

# get file directory 
model_directory = find_directory('models', __file__)
data_directory = find_directory('Data', __file__)

# initialize model
E0, S0, all_components = BasicModel(model_directory=model_directory, flux_weighting=True).load_ES_models()

# FBA alpha - m1
df = get_div_obj_df([E0, S0], target_obj_val=0.5, potential_genes=['folA'], precision=2)
df.to_csv(os.path.join(data_directory, 'alpha_table_m1.csv'))

# TODO: refine nonessential entry and m3 table