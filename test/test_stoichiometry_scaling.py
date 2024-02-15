import logging
import pandas as pd

from scarcc.preparation.perturbation import alter_Sij, get_alphas_from_tab
# from scarcc.preparation.perturbation import iter_species
from scarcc.preparation.metabolic_model import BasicModel

logging.basicConfig(filename='test_stoi_scaling.log', 
                    level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')

E0, S0, all_components = BasicModel(flux_weighting=True).load_ES_models()
alpha_table = pd.read_csv('../Data/alpha_table_m1.csv', index_col=0)

def iter_species(models,f,*args,**kwargs): 
    r_object = list()
    if type(models) is not zip:
        for model in models:
            print(*args,**kwargs)

            r_object.append(f(model,*args,**kwargs))
    else:
        for model, *extra_objects in models: 
            r_object.append(f(model,extra_objects,*args,**kwargs))
    return(r_object) 

def change_model_stoichiometry(model, alphas):
    with model:
        changed_rxns = alter_Sij(model, alphas=alphas, genes='folA', ko=False)
        print(changed_rxns)
        logging.debug('For model %s with alpha %s, Scaled reactions: %s', model.id, alphas[0],'\n'.join([str(model.reactions.get_by_id(rxn)) for rxn in changed_rxns]))
    
species_list = [E0, S0]
alpha_list = [get_alphas_from_tab(model, ['folA'], alpha_table) for model in species_list]
print(alpha_list)
zipped_sp_alpha = zip([E0, S0], alpha_list)
iter_species(zipped_sp_alpha, change_model_stoichiometry)