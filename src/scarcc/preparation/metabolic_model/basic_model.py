"""Basic model class for loading cross-feeding species with customized model adjustment"""

from dataclasses import dataclass
from typing import Any, Dict
import itertools
import logging
import os

import cobra


# from .core.medium import initialize_medium, change_medium
# from .core.component import get_all_components
# from .core.flux_weighting import weight_carbon_byproduct, rename_ac_in_bulk

from .core import (initialize_medium, change_medium,
                   get_all_components,
                   weight_carbon_byproduct, rename_ac_in_bulk)

logger = logging.getLogger(__name__)

@dataclass
class BasicModel:
    """Class for loading cross-feeding species with customized model adjustment
    
    Args:
        E0: cobra.Model
        S0: cobra.Model
        all_components: dict, Storing all components(Gene, Metabolite, Reaction) associated with metabolic networks
        E_carbon_byproduct_reuptake: Boolean, default True, allow for reuptake of carbon byproduct in E0
                                        pFBA galactose reuptake adds noise in growth rate estimation
        flux_weighting: Boolean, default None, flux weighting lower relative importance in cost of acetate secretion in pFBA
                                                                increase relative importance in cost of galactose secretion in pFBA
        ac_scale: float, default 10, flux weight 1/10
        gal_scale: float, default 3, flux weight 3
        c_limiting_conc: float, default 2, carbon limiting as proxy for coculture uptake rate
        met_limiting_conc: float, default 10, metabolite unlimited for E0        

    Returns:
        _type_: _description_
    """
    
    E0: cobra.Model = None
    S0: cobra.Model = None
    model_directory: str = '../models'
    all_components: Dict[str, Any] = None
    E_carbon_byproduct_reuptake: bool = True
    flux_weighting: bool = None
    ac_scale: float = 10
    gal_scale: float = 3
    c_limiting_conc = 2 # corresponding to maximum uptake rate in carbon limited
    met_limiting_conc = 10 # DO not search monoculture alpha with met_limiting_conc
    
    def __post_init__(self):
        pass

    def set_ES_medium(self):
        """Set the medium for E0 and S0 models"""
        nutrient_medium = initialize_medium()
        self.E0.medium, self.S0.medium = nutrient_medium, nutrient_medium        
        change_medium(self.E0, ['EX_lcts_e', 'EX_met__L_e'], [self.c_limiting_conc, self.met_limiting_conc])
        change_medium(self.S0, 'EX_gal_e', self.c_limiting_conc)
        return None

    def implement_flux_weighting(self):
        """Flux weighting for carbon byproduct secretion in E0
            Reduce flux count during minimization of fluxes in pFBA for acetate
            Increase flux count during minimization of fluxes in pFBA for galactose
            renaming of ac_p and ac_e to bulk_ac_p and bulk_ac_e in both models to allow for valid exchange
                
                Note: COMETS implementation do not allow for weighting exchange flux directly(EX_ac_e: 10*ac -> ),
                To reduce flux quantity the same way as weighting exchange flux, I equate 10 ac_p and 1 bulk_ac_p, 10 ac_e and 1 bulk_ac_e
                Producing equivalent exchange fluxes weighting, S0 model adjustment for applying the same reaction and metabolite change
        """
        add_metab_quantity = 1-1/self.ac_scale
        self.E0.reactions.ACt2rpp.add_metabolites({self.E0.metabolites.ac_p: add_metab_quantity}) # 0.1 ac_p <-> ac_c
        self.S0.reactions.ACt2rpp.add_metabolites({self.S0.metabolites.ac_p: add_metab_quantity})
        weight_carbon_byproduct(self.E0, self.S0, self.all_components, ac_scale=self.ac_scale, gal_scale=self.gal_scale)

        for model, metabolite in itertools.product([self.E0, self.S0], ['ac_p', 'ac_e']):
            rename_ac_in_bulk(model, metabolite)
            
    def load_ES_models(self):
        """Load E0 and S0 models with customized """
        if self.E0 is None:
            self.E0 = cobra.io.read_sbml_model(os.path.join(self.model_directory, "iML1515_E0.xml"))
            self.S0 = cobra.io.read_sbml_model(os.path.join(self.model_directory, "STM_v1_0_S0.xml"))
            self.E0.id, self.S0.id = 'E0', 'S0'
        if self.E0 is None or self.S0 is None:
            logging.warning('%s and %s models are not loaded', self.E0.id, self.S0.id)
        else:
            logging.debug('%s and %s models are loaded', self.E0.id, self.S0.id)
        # galactose reuptake in Diauxic environment result in disturbance in growth rate estimate after lcts being used up
        if self.E_carbon_byproduct_reuptake is False:
            logging.debug('GAL reuptake is turned off')
            self.E0.reactions.GALtex.upper_bound = 0
            self.E0.reactions.ACtex.upper_bound = 0
        self.set_ES_medium()
        self.all_components = get_all_components(self.E0, self.S0)
        if self.flux_weighting is True:
            logging.debug('Flux weighting is turned on')
            self.implement_flux_weighting()
        return self.E0, self.S0, self.all_components
    