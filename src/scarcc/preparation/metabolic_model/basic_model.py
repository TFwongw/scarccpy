import cobra
from cobra.io import load_model
import cometspy as c
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from functools import partial
import re
from dataclasses import dataclass, field
from typing import Callable, List, Any, Dict
import sys
import itertools
import json
import logging
import multiprocessing
from functools import partial
import copy
import re

import dataclasses

# from ..util.util import rename_columns
from .core.component import get_links_component, get_component, get_all_components
from .core.medium import initialize_medium, change_medium

from .core.flux_weighting import weight_carbon_byproduct, rename_ac_in_bulk

import logging
logger = logging.getLogger(__name__)


@dataclass
class basic_model:
    E0: cobra.Model = None
    S0: cobra.Model = None
    all_components: Dict[str, Any] = None
    Egal_reuptake: bool = True
    flux_weighting: bool = None
    ac_scale: float = None
    gal_scale: float = None

    def __post_init__(self):
        pass

    def set_ES_medium(self):
        nutrient_medium = initialize_medium()
        self.E0.medium, self.S0.medium = nutrient_medium, nutrient_medium
        c_limiting_conc = 10 # corresponding to maximum uptake rate in carbon limited
        met_limiting_conc = 10 # DO not search momoculture alpha with met_limitng_conc
        
        change_medium(self.E0, ['EX_lcts_e', 'EX_met__L_e'], [c_limiting_conc , met_limiting_conc])
        change_medium(self.S0, 'EX_gal_e', c_limiting_conc)
        return None

    def implement_flux_weighting(self):
        weight_carbon_byproduct(self.E0, self.S0, self.all_components, ac_scale=self.ac_scale, gal_scale=self.gal_scale)
        self.E0.reactions.ACt2rpp.add_metabolites({self.E0.metabolites.ac_p: .9}) # 0.1 ac_p <-> ac_c
        self.S0.reactions.ACt2rpp.add_metabolites({self.S0.metabolites.ac_p: .9})

        [rename_ac_in_bulk(model, metabolite) for model, metabolite in itertools.product([self.E0, self.S0], ['ac_p', 'ac_e'])] # rename ac to bulk_ac for both metab in both models
        return None

    def load_ES_models(self):
        if self.E0 is None:
            self.E0 = cobra.io.read_sbml_model("../models/iML1515_E0.xml")
            self.S0 = cobra.io.read_sbml_model("../models/STM_v1_0_S0.xml")
            self.E0.id, self.S0.id = 'E0', 'S0' 
        
        logging.debug(f'{self.E0.id} and {self.S0.id} models are loaded')
        
        # galactose reuptake in Diauxic environment result in disturbance in growth rate estimate after lcts being used up
        if self.Egal_reuptake is False:
            logging.debug('GAL reuptake is turned off')
            self.E0.reactions.GALtex.upper_bound = 0
            self.E0.reactions.ACtex.upper_bound = 0
        self.set_ES_medium()
        self.all_components = get_all_components(self.E0, self.S0)
        if self.flux_weighting is True:
            logging.debug('Flux weighting is turned on')
            self.implement_flux_weighting()
        
        if self.E0 is None or self.S0 is None:
            logging.debug(f'{self.E0.id} and {self.S0.id} models are not loaded')
        else: logging.debug(f'{self.E0.id} and {self.S0.id} models second')
        return self.E0, self.S0, self.all_components

# def rename_ac_in_bulk(model, id):
#         new_id = f'bulk_{id}'
#         model.metabolites.get_by_id(id).id = new_id
#         if 'ac' in id:
#             model.metabolites.get_by_id(new_id).name = 'bulk Acetate'
#             model.metabolites.get_by_id(new_id).formula = 'C20H30O20'
#         if '_e' in id:
#             model.reactions.get_by_id(f'EX_{id}').id = f'EX_{new_id}'

# def weight_carbon_byproduct(E0, S0, all_components, ac_scale=None, gal_scale=None):
#     if gal_scale is not None:
#         query_gal = ['gal_p','gal_c']
#         if gal_scale > 1:
#             gal_scale = -1*(1-1/gal_scale) 
#         for rxn in get_links_component(E0, query_gal, all_components, id_only=False, is_prod_only=True): # ['GALt2pp', 'GAL1PPpp', 'GALabcpp', 'GALS3', 'LACZpp', 'GALM2pp', 'LACZ'] are scaled
#             # TODO: only scale LACZpp is desirable, make gal secretion similar to knockout
#             metab_to_scale = rxn.metabolites
#             rxn.add_metabolites({k:v*gal_scale for k,v in metab_to_scale.items()})
#     # E0.reactions.GALtex.knock_out()
    
#     if ac_scale is not None: # flux do not count in minflux
#         add_scale = ac_scale-1 if ac_scale>=1 else ac_scale
        
#         # 0.1 ac_p + 10 h_p <=> 10 ac_c + 10 h_c
#         for rxn in [ 'ACt2rpp']: 
#             if isinstance(rxn, str):
#                 rxn = get_component(E0, rxn, all_components)
#             if not 'EX_' in rxn.id:
#                 metab_to_scale = rxn.metabolites
#                 rxn.add_metabolites({k:v*add_scale for k,v in metab_to_scale.items()}) 

# ##        
# # from .iter_species import iter_species
# from ..perturb.iter_species import iter_species
# from ..perturb.stoichiometry_scaling import get_alphas_from_tab

# from ..sim_engine.simulation_builder import create_c
# # from ..sim_engine.result_processing import 

# from ..sim_engine.simulation_builder import get_BM_df

# def remove_Zero_col(df): # extend N differ than 0 
#     return(df.loc[:, ((df !=0) & (df.notnull())).any(axis=0)]) # ignore NA entry 

# # %%
