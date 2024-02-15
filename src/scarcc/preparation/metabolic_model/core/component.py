"""This module contains functions for retriving components of a model"""
import re
import cobra
from typing import Union, List, Dict
# from scarcc.util.util import convert_arg_to_list

def get_gene_id(model: "cobra.Model", gene_name: str) -> str:
    """Query the id of the gene
    
    Parameters
    ----------
    model: cobra.Model
    gene_name: str
    """
    for i in model.genes:
        if(i.name == gene_name):
            return(i.id)

def get_all_components(E0, S0) -> Dict[str, Dict[str, List]]: 
    """Returns all the components of the model
    
    Parameters
    ----------
    E0: cobra.Model
    S0: cobra.Model

    Returns
    -------
     Structured Dictionary:
        all_components = {'metabolites': {'E0': [met1, met2, ...], 'S0': [met1, met2, ...]},
                        'genes': {'E0': [gene1, gene2, ...], 'S0': [gene1, gene2, ...]}, 
                        'reactions': {'E0': [rxn1, rxn2, ...], 'S0': [rxn1, rxn2, ...}}
        }
    """

    all_metabolites = {
        model.id: model.metabolites for model in [E0, S0]
    }
    all_genes = {
        model.id: 
        [gene.name for gene in model.genes if type(gene) != str]
        for model in [E0, S0]
    }

    all_reactions = {
        model.id: model.reactions for model in [E0, S0]
    }

    all_components = {'metabolites': all_metabolites,
                    'genes': all_genes, 
                    'reactions': all_reactions}
    return all_components

def substrate_only(rxn: "Reaction", met: Union["Metabolite" , str]) -> bool:
    """Returns True if the metabolite is a substrate but not product in the reaction
    
    Parameters
    ----------
    rxn: cobra.Reaction
    met: cobra.Metabolite
    
    Returns
    -------
    bool
    """
    return rxn.lower_bound>=0 and met in rxn.reactants        

def substrate_any(rxn: "Reaction", met: Union["Metabolite" , str]):
    """Returns True if the metabolite can be a substrate in the reaction"""
    return (rxn.lower_bound>=0 and met in rxn.reactants) or rxn.lower_bound<0

def product_only(rxn: "Reaction", met: Union["Metabolite" , str]):
    """Returns True if the metabolite is a product but not substrate in the reaction
       Excludes secretion fluxes"""
    if re.search(r"_e$", rxn.id): # secretion flux
        return True
    return (rxn.lower_bound>=0 and met in rxn.products) or (rxn.upper_bound<=0 and met in rxn.reactants)   
    
def get_component(model, query, all_components):
    """Returns the component of the model that matches the query"""
    ID = model.id.split('.')[0]
    if query in all_components['genes'][ID]:
        return model.genes.get_by_id(get_gene_id(model, query))
    if query in all_components['metabolites'][ID]:
        return model.metabolites.get_by_id(query)
    return model.reactions.get_by_id(query)
    
def get_links_component(model: "Model", query: Union["Metabolite" , "Reaction", "Gene", str], all_components: Dict, id_only: bool = True, 
                        is_sub_only: bool = None, is_any_sub: bool = None, is_any_prod: bool = None, is_prod_only: bool = None):
    """Returns the components linked to the query
    
    Parameters
    ----------
    model: cobra.Model
        The model to search

    query: str or list
        The component to search for, can be a list of genes, metabolites or reactions
    
    all_components: dict
        The dictionary containing all the components of the model
    
    id_only: bool
        If True, returns the id of the components
    
    is_sub_only: bool
        If True, returns only the reactions where the query is a substrate
    
    is_any_sub: bool
        If True, returns the reactions where the query is a substrate
    
    is_any_prod: bool
        If True, returns the reactions where the query is a product
    
    is_prod_only: bool
        If True, returns only the reactions where the query is a product
        
    Returns
    -------
    list
        A list of the components linked to the query
        [component1, component2, ...] or [component1.id, component2.id, ...]
    """

    def get_result(model, query, all_components):        
        # query for single component
        query = get_component(model, query, all_components)
        if isinstance(query, cobra.Reaction):
            result = set(query.metabolites)
        else:
            result = query.reactions
            if isinstance(query, cobra.Metabolite):
                if is_sub_only is not None:
                    result = set([rxn for rxn in result if substrate_only(rxn, query)==is_sub_only])
                elif is_any_sub is not None or is_any_prod is not None:
                    result = set([rxn for rxn in result if substrate_any(rxn, query)==is_any_sub])
                elif is_prod_only is not None:
                    result = set([rxn for rxn in result if product_only(rxn, query)==is_prod_only])
        return result
    
    result=set()
    for ele in convert_arg_to_list(query):
        result.update(get_result(model, ele, all_components))
    if id_only:
            result = set([ele.id for ele in result])
    return list(result)
