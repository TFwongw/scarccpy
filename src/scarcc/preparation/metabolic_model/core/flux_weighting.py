"""Flux weighting as an approach to refine carbon byproduct secretion flux for E.coli in pfba"""

import logging
logger = logging.getLogger(__name__)

from . import get_links_component, get_component

def rename_ac_in_bulk(model, id):
    """Change name of acetate metabolite to a single bulk acetate
    
    Parameters
    ----------
    model: cobra.Model
    id: str
        id of metabolite to be renamed
    
    Returns
    -------
    None
    """
    
    new_id = f'bulk_{id}'
    model.metabolites.get_by_id(id).id = new_id
    if 'ac' in id:
        model.metabolites.get_by_id(new_id).name = 'bulk Acetate'
        model.metabolites.get_by_id(new_id).formula = 'C20H30O20'
    if '_e' in id:
        model.reactions.get_by_id(f'EX_{id}').id = f'EX_{new_id}'

def weight_carbon_byproduct(E0, S0, all_components, ac_scale=None, gal_scale=None):
    """Weight carbon byproduct secretion flux for E.coli with specific scales
    If a scale is None, no change is made
    Not changing functionality of metabolic model unlike stoichiometric scaling
    
    Parameters
    ----------
    E0: cobra.Model
    S0: cobra.Model
    all_components: dict
    ac_scale: float
        scale need to be larger than 1, reduce in flux accounted in pfba optimization procedure
    gal_scale: float
        scale need to be less than 1, increase in flux accounted in pfba optimization procedure

    Returns
    -------
    None
    """
    if gal_scale is not None:
        gal_scale = -1*(1-1/gal_scale) if gal_scale > 1 else -1*(1-gal_scale) # ensure increase in cose
        
        for rxn in [E0.reactions.get_by_id('LACZpp')]:
            # TODO: only scale LACZpp is desirable, make gal secretion similar to knockout
            metab_to_scale = rxn.metabolites
            rxn.add_metabolites({k:v*gal_scale for k,v in metab_to_scale.items()}) # similar toE0.reactions.GALtex.knock_out()
    
    if ac_scale is not None: # flux do not count in minflux
        add_scale = ac_scale-1 if ac_scale>=1 else ac_scale
        
        # ac_p + 10 h_p <=> 10 ac_c + 10 h_c
        for rxn in ['ACt2rpp']:
            if isinstance(rxn, str):
                rxn = get_component(E0, rxn, all_components)
            if not 'EX_' in rxn.id:
                metab_to_scale = rxn.metabolites
                rxn.add_metabolites({k:v*add_scale for k,v in metab_to_scale.items()})
