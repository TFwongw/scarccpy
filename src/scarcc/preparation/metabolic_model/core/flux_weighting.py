import logging
logger = logging.getLogger(__name__)

from . import get_links_component, get_component

def rename_ac_in_bulk(model, id):
    
    new_id = f'bulk_{id}'
    model.metabolites.get_by_id(id).id = new_id
    if 'ac' in id:
        model.metabolites.get_by_id(new_id).name = 'bulk Acetate'
        model.metabolites.get_by_id(new_id).formula = 'C20H30O20'
    if '_e' in id:
        model.reactions.get_by_id(f'EX_{id}').id = f'EX_{new_id}'

def weight_carbon_byproduct(E0, S0, all_components, ac_scale=None, gal_scale=None):
    if gal_scale is not None:
        query_gal = ['gal_p','gal_c']
        if gal_scale > 1:
            gal_scale = -1*(1-1/gal_scale) 
        for i, rxn in enumerate(get_links_component(E0, query_gal, all_components, id_only=False, is_prod_only=True)): # ['GALt2pp', 'GAL1PPpp', 'GALabcpp', 'GALS3', 'LACZpp', 'GALM2pp', 'LACZ'] are scaled
            # TODO: only scale LACZpp is desirable, make gal secretion similar to knockout
            metab_to_scale = rxn.metabolites
            rxn.add_metabolites({k:v*gal_scale for k,v in metab_to_scale.items()})
            logger.debug('%s th reaction in gal flux weighted: \n %s', i, str(rxn))
    # E0.reactions.GALtex.knock_out()
    
    if ac_scale is not None: # flux do not count in minflux
        add_scale = ac_scale-1 if ac_scale>=1 else ac_scale
        
        # 0.1 ac_p + 10 h_p <=> 10 ac_c + 10 h_c
        for rxn in ['ACt2rpp']: 
            if isinstance(rxn, str):
                rxn = get_component(E0, rxn, all_components)
            if not 'EX_' in rxn.id:
                metab_to_scale = rxn.metabolites
                rxn.add_metabolites({k:v*add_scale for k,v in metab_to_scale.items()}) 
