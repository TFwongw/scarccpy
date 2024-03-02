"""functions for stoichiometry scaling"""

import pandas as pd
from typing import List
from scarcc.utils import convert_arg_to_list
from scarcc.preparation.metabolic_model import get_gene_id

def scale_reaction(model: "Model", reaction_id: str, alpha: float, direction='forward'):
    """
    scale a reaction by alpha, change the direction of the reaction if needed
    
    Parameters
    ----------
    
    model: cobra.Model
    reaction_id: str
    alpha: float
    direction: str
        'forward': scale product in the reaction 
        'backward': scale reactant in the reaction
    """
    if direction == 'forward':
        model.reactions.get_by_id(reaction_id).lower_bound = 0
        dict_product = dict((k, v) for k, v in model.reactions.get_by_id(reaction_id).metabolites.items() if v >= 0)    # obtain 'end product'
    else:  # reverse reaction
        model.reactions.get_by_id(reaction_id).upper_bound = 0
        dict_product = dict((k, v) for k, v in model.reactions.get_by_id(reaction_id).metabolites.items() if v <= 0)
    if isinstance(alpha, pd.Series):
        alpha = float(alpha)

    for product, unit in dict_product.items():  # scale corresponding metabolite units involved in the reaction
        model.reactions.get_by_id(reaction_id).add_metabolites({product: -unit*(1-1/alpha)}) # only change unit of unidirection metabolite

def separate_reaction(model: "Model", reaction_id: str, alpha: float) -> List:
    """
    decompose reversible reaction into forward and backward
    scale each end product by 1/alpha

    Parameters
    ----------
    model: cobra.Model
    reaction_id: str
    alpha: float
    
    Returns
    -------
    list of reaction ids that are regulated by the target genes
    """
    (lb, ub) = model.reactions.get_by_id(reaction_id).bounds 
    rxn_ids = [reaction_id] #? 
    if(lb < 0 and ub !=0): # only perform if reaction is bidirectional
        rev_reaction = model.reactions.get_by_id(reaction_id).copy() # copy of target reaction 
        rev_reaction.id = f'{reaction_id}_v1' # redefine id for copy of reaction
        model.add_reactions([rev_reaction]) # add to model
    
        scale_reaction(model, rev_reaction.id, alpha, direction='backward')      
        
        rxn_ids.append(rev_reaction.id) # list of id for reaction and reaction_v1 
    scale_reaction(model, reaction_id, alpha, direction='forward')
    return(rxn_ids)

def alter_Sij(model: "Model", alphas: float = 1, genes: str = 'folA', ko=False):
    """Change stoichiometry of reactions encoded by target genes to alpha * original stoichiometry
    
    Parameters
    ----------
    model: cobra.Model
    alphas: float or list of float
        scaling factor for target genes
    genes: str or list of str
        gene id or gene name
    ko: bool
        knock out the reaction if True, otherwise scale stoichiometry
    
    Returns
    -------
    list of reaction ids that are regulated by the target genes
    """
    # get objective value for corresponding alpha
    alphas= convert_arg_to_list(alphas[0]) if isinstance(alphas, list) and len(alphas)==1 else convert_arg_to_list(alphas)  # unlist one layer from zip comprehension 
    genes = convert_arg_to_list(genes)
    
    genes_dict = {gene: alpha for gene, alpha in zip(genes, alphas)}
    genes_sorted = sorted(genes_dict.items(), key=lambda x:x[1], reverse=True) #sort by magnitude of alpha
    rxn_ids = list() # store list of id of reaction and reaction_v1 that regulated by the same gene 
    for current_gene, alpha in genes_sorted:
        current_gene = current_gene.split('_')[0] # for step_alpha
        for rxn in model.genes.get_by_id(get_gene_id(model, current_gene)).reactions:
            if ko:
                model.reactions.get_by_id(rxn.id).knock_out()
            elif (rxn.id not in rxn_ids):
                rxn_ids.extend(separate_reaction(model, rxn.id, alpha))# copy of reaction, forward_terminal_change = True
    return(rxn_ids)

def get_alphas_from_tab(model: "Model", genes: List, alpha_table: pd.DataFrame) -> List: # the returned df of one gcomb index only pass as series
    """Get alpha values from alpha_table for target genes

    Parameters
    ----------
    model: cobra.Model
    genes: list
        gene id or gene name
        alpha_table: pd.DataFrame
            alpha values for genes and reactions
    
    Returns
    -------
    list of alpha values for target genes
    """
    if isinstance(alpha_table, pd.Series):
        alpha_table = alpha_table.to_frame().T
    genes = convert_arg_to_list(genes)
    alphas = [alpha_table.loc[gene ,model.id] for gene in genes]
    return alphas