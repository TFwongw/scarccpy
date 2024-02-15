

def scale_reaction(model, reaction_id, alpha, direction='forward'):
    if direction == 'forward':
        model.reactions.get_by_id(reaction_id).lower_bound = 0
        dict_product = dict((k, v) for k, v in model.reactions.get_by_id(reaction_id).metabolites.items() if v >= 0)    # obetain 'end product'
    else:  # reverse reaction
        model.reactions.get_by_id(reaction_id).upper_bound = 0
        dict_product = dict((k, v) for k, v in model.reactions.get_by_id(reaction_id).metabolites.items() if v <= 0)
    if isinstance(alpha, pd.Series):
        alpha = float(alpha)

    for product, unit in dict_product.items():  # scale corresponding metabolite units involved in the reaction
        model.reactions.get_by_id(reaction_id).add_metabolites({product: -unit*(1-1/alpha)}) # only change unit of unidirection metabolite
    return None

def separate_reaction(model, reaction_id, alpha):
    """
    decompose reversible reaction into forward and backward
    scale each end product by 1/alpha
    """
    (lb, ub) = model.reactions.get_by_id(reaction_id).bounds 
    rct_ids = [reaction_id] #? 
    if(lb < 0 and ub !=0): # only perform if reaction is bidirectional
        rev_reaction = model.reactions.get_by_id(reaction_id).copy() # copy of target reaction 
        rev_reaction.id = f'{reaction_id}_v1' # redefine id for copy of reaction
        model.add_reactions([rev_reaction]) # add to model
    
        scale_reaction(model, rev_reaction.id, alpha, direction='backward')      
        
        rct_ids.append(rev_reaction.id) # list of id for reaction and reaction_v1 
    scale_reaction(model, reaction_id, alpha, direction='forward')
    return(rct_ids)

def alter_Sij(model, alphas = 1, genes = 'folA', ko=False):  
    # get objective value for corresponding alpha
    alphas= convert_arg_to_list(alphas[0]) if type(alphas) is list and len(alphas)==1 else convert_arg_to_list(alphas)  # unlist one layer from zip comprehension 
    genes = convert_arg_to_list(genes)
    
    genes_dict = {gene: alpha for gene, alpha in zip(genes, alphas)}
    genes_sorted = sorted(genes_dict.items(), key=lambda x:x[1], reverse=True) #sort by magnitude of alpha
    rct_ids = list() # store list of id of reaction and reaction_v1 that regulated by the same gene 
    for current_gene, alpha in genes_sorted:
        current_gene = current_gene.split('_')[0] # for step_alpha
        for rct in model.genes.get_by_id(get_gene_id(model, current_gene)).reactions: 
            if ko:
                model.reactions.get_by_id(rct.id).knock_out()
            elif (rct.id not in rct_ids):
                rct_ids.extend(separate_reaction(model, rct.id, alpha))# copy of reaction, forward_terminal_change = True
    return(rct_ids) 




def get_alphas_from_tab(model, genes: list, alpha_table): # the returned df of one gcomb index only pass as series
    if isinstance(alpha_table, pd.Series):
        alpha_table = alpha_table.to_frame().T
    genes = convert_arg_to_list(genes)
    alphas = [alpha_table.loc[gene ,model.id] for gene in genes]
    return alphas 

