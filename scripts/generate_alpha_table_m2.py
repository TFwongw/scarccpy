# This script generates the alpha table using alpha_table_m1 and 
# keeping species alpha ratio indicated there to search for the alpha pairs 
# that corresponding to IC50 values in coculture growth. 

import os
import pandas as pd

from scarcc.preparation.metabolic_model import BasicModel
from scarcc.preparation.find_directory import find_directory

# read potential genes
potential_genes = ['folP','folA']

# get file directory 
model_directory = find_directory('models', __file__)
data_directory = find_directory('Data', __file__)

# initialize model
E0, S0, all_components = BasicModel(model_directory=model_directory, flux_weighting=True).load_ES_models()

# Method 2 coculture alpha search
from scarcc.preparation.perturbation.alpha_finder.coculture import CocultureAlphaFinder, run_coculture_search_mp
import cometspy as c

alpha_table = pd.read_csv(os.path.join(data_directory, 'alpha_table_m1.csv'), index_col=0)
initial_pop = 1e-8

p = c.params()
p.set_param("defaultKm", 0.00001) # M 
p.set_param("defaultVmax", 10) #mmol/gDw/hr
p.set_param("maxCycles", 180)
p.set_param("timeStep", 1) 
p.set_param('writeFluxLog', True)
p.set_param('writeMediaLog', True)

obj_style = 'MAX_OBJECTIVE_MIN_TOTAL'

def get_gr_Normal(): # Growth rate estimation without drug
    AF = CocultureAlphaFinder(model=[E0, S0], search_alpha=None, current_gene = 'Normal', p=p,
                            exp_leap= 1.3, alpha_table=alpha_table, carbon_source_val=.1,
                            target_obj_val=.5, precision=1, initial_pop=initial_pop, obj_style=obj_style)
    gr_Normal = AF.calculate_gr_Normal()
    print(f'---------GRNROMAL {gr_Normal}')
    return gr_Normal
    
gr_Normal = get_gr_Normal()
search_alpha = None

kwargs = {
    'exp_leap': 1.5, # for coculture gr inherently > 1, require large leap 
    'alpha_table': alpha_table,
    'p': p, 
    'carbon_source_val': .1,
    'target_obj_val': .5,
    'precision': 2, 
    'gr_Normal': gr_Normal,
    'initial_pop': initial_pop, 
    'obj_style': obj_style
}

# alpha_table_m2 get stored in Data directory
opt_alpha_list, trace_dict = run_coculture_search_mp(potential_genes, data_path=data_directory,file_suffix='m2', n_processor=10,
                                                                model=[E0, S0], search_alpha=search_alpha,
                                                                **kwargs)
