# Find alpha from 

from dataclasses import dataclass, field
import logging
import pandas as pd
from typing import List

from scarcc.preparation.perturbation import alter_Sij, iter_species
from scarcc.util import convert_arg_to_list
from .alpha_finder import AlphaFinderConfig

logger = logging.getLogger(__name__)

potential_genes = [
    'glyA', 'gltA', 'tktA', 'dapF', 'dapB', 'acnB', 'pgk', 'talB', 'mrcA', 'pyrE', 'dapD', 'pfkA', 'gdhA', 'folA', 'mrdA', 'thrB',
    'dapA', 'serC', 'argD', 'thrC', 'aceF', 'pykF', 'dadX', 'folC', 'pyrD', 'trpA', 'serB', 'fbp', 'eno', 'pgi', 'pheA',
    'gcvH', 'gnd', 'murA', 'aroA', 'guaB', 'glnA', 'yrbG', 'folP', 'purU', 'serA', 'gltD', 'purT', 'ackA', 'purN', 'rffG',
    'gapA'
]

@dataclass(kw_only=True)
class MonocultureAlphaFinder(AlphaFinderConfig):
    model : 'cobra.Model' # single cobra model
    exp_leap : int = 2
    trace_obj_val : List = field(default_factory=list)
    norm_obj_val : float = None # gr_Nomral in coculture case
    response_record = None # for iteratice uptate in checkerboard
    opt_df : pd.DataFrame = None

    def __post_init__(self):
        # super().__post_init__()

        self.is_monoculture = True
        if not self.acceptance_threshold_lower:
            self.acceptance_threshold_lower = .85
        if not self.acceptance_threshold_upper:
            self.acceptance_threshold_upper = 1.1
        self.norm_obj_val = self.model.slim_optimize()
        if not self.response_record:       
            self.response_record = {self.current_gene: {self.model: {'ICX': convert_arg_to_list(self.target_obj_val), 'response': {}}}}

    def fill_response_record(self, obj_val, is_lowest_abs_diff=False):
        obj_val_interval = float(format(obj_val, '.2f'))
        record = self.response_record[self.current_gene][self.model]['response'].get(obj_val_interval) # closest alpha in the obj_interval
        interval_abs_diff = abs(obj_val - obj_val_interval) 

        if record and (record['interval_abs_diff'] >= interval_abs_diff): # erase record, update to current alpha
            record = None  
        if record and (obj_val_interval > .99) and (record['search_alpha']>self.search_alpha): # IC0 update most large alpha
            record = None

        if self.opt_df is not None:
            previous_min = min([ele['target_abs_diff'] for ele in 
                    self.response_record[self.current_gene][self.model]['response'].values()])
            if interval_abs_diff <= previous_min:
                is_lowest_abs_diff = True

        if not record:
            self.response_record[self.current_gene][self.model]['response'][obj_val_interval] = {
                'search_alpha' : self.search_alpha,
                'precise_obj_val' : obj_val, 
                'interval_abs_diff' : interval_abs_diff,
                'target_abs_diff' : abs(obj_val - self.target_obj_val)}
        return is_lowest_abs_diff

    def eval_alpha_fun(self):
        # print(self.alpha_ub)
        _, obj_val, summary_df = Sij_biomass(self.model, self.search_alpha, self.current_gene)
        
        if(self.opt_df is None and summary_df['Net_Flux'][0]<1e-8): # reinitialize if the first alpha to search gives zero flux
            self.search_alpha = 1+2e-3
            _, obj_val, summary_df = Sij_biomass(self.model, self.search_alpha, self.current_gene)
        self.trace_obj_val.append(obj_val)
        obj_val = obj_val/ self.norm_obj_val # to normalized 
        self.is_new_ub = (summary_df['Net_Flux'][0]<1e-8 or summary_df['Net_Flux'][0]>100 or 
                        obj_val<self.target_obj_val*self.acceptance_threshold_lower or obj_val>1.2) # overinhibition, search alpha too high 
        
        self.found_alpha = self.is_found(self.search_alpha, self.alpha_lb, self.alpha_ub, self.precision)
        
        # update optimal df
        net_flux_req = (summary_df['Net_Flux'][0]>0 or obj_val>self.norm_obj_val*.1)
        is_lowest_abs_diff = self.fill_response_record(obj_val)

        obj_req = is_lowest_abs_diff
        if (net_flux_req and obj_req) or (self.opt_df is None): # store only qualified alpha
            self.opt_df = summary_df
        
    def out_fun(self):
        opt_df = self.opt_df
        opt_df['is_growth_switch'] = self.classify_growth_switch()
        logger.debug('Gene: %s with alpha %s and obj_val %s', self.current_gene, opt_df['div_opt_alpha'][0], opt_df[f'div_opt_obj_val'][0])
        opt_df.insert(loc=2, column='Percent_target_obj_val', value=opt_df['div_opt_obj_val']/self.norm_obj_val)
        opt_df.columns = [f'{self.model.id}_'+element for element in opt_df.columns]
        self.opt_df = opt_df
        return (opt_df[f'{self.model.id}_div_opt_alpha'],opt_df[f'{self.model.id}_div_opt_obj_val'],opt_df) # alpha_feasible, and upper bound of alpha

def get_summary_df(model, alpha, obj_val, rct_ids = 'DHFR', sol_sol = None): 
#     summarize solutions from optimization and alpha used  
# expect direction opposite with zero flux, to be consistent with FVA bound 
    if type(alpha) != int and type(alpha) != float:
        alpha = str(alpha)
    if sol_sol is not None:
        rct_dir, rct_rev, flux_dir, flux_rev = list(), list(), list(), list()
        for current_rct_id in rct_ids:
            append_str = str(model.reactions.get_by_id(current_rct_id))
            if ('_v1' not in current_rct_id): # direction align with FVA
                rct_dir.extend([append_str])
                try:
                    flux_dir.extend([round(model.reactions.get_by_id(current_rct_id).flux,5)])
                except: 
                    flux_dir.extend([sol_sol[current_rct_id]])
            else: # direction opposite
                rct_rev.extend([append_str])        
                try:
                    flux_rev.extend([round(model.reactions.get_by_id(current_rct_id).flux,5)])
                except:
                    flux_rev.extend([sol_sol[current_rct_id]])
    #     rct_rev = fix_string(rct_rev)
    #     rct_dir = fix_string(rct_dir)
        net_flux = sum(abs(element) for element in set(flux_dir) | set(flux_rev))
        net_flux_I = 'Zero Flux' if net_flux ==0 else 'Net Flux'

    #     if type(alpha) is list or type(alpha) is tuple:

        summary_df = pd.DataFrame({f'div_opt_alpha': alpha,
                                    f'div_opt_obj_val': obj_val,
                                    f'FVAdir_Reactions_id': ', '.join(rct_ids),
                                    f'FVAdir_Reactions_involved': ' '.join(rct_dir),                                    
                                    f'FVAdir_Flux_values':[flux_dir],
                                    f'FVAopposite_Reactions_involved': ', '.join(rct_rev),     
                                    f'FVAopposite_Flux_values':[flux_rev],
                                    f'Net_Flux_Boolean':[net_flux_I],
                                    f'Net_Flux':[net_flux]}
                                    )
#     else:
#         summary_df = pd.DataFrame([{'div_opt_obj_val': obj_val}])        
    return(summary_df)


def Sij_biomass(model, alphas = 1, genes = 'folA', slim_opt = False): 
    with model:
        rct_ids = alter_Sij(model, alphas, genes) # same

        if (slim_opt == False):
            sol_sol = model.optimize()
            obj_val = sol_sol.objective_value 
            summary_df = get_summary_df(model, alphas, obj_val, rct_ids, sol_sol.fluxes)  
            return(alphas, obj_val, summary_df) 
        else:
#             return(pd.DataFrame([{f'div_opt_obj_val': model.slim_optimize()}]))
            return(model.slim_optimize())

def get_single_div_obj_df(model, target_obj_val, first_n_gene=None, alpha_df = None, search_alpha = None,precision=6,
                        potential_genes=potential_genes, acceptance_threshold_upper=1.1, acceptance_threshold_lower=.95): # media not required for GR

#     # generate full biomass dataframe for a single species
#     # given alpha or automatic search for optimal with the given objective value
#     if not first_n_gene: first_n_gene = len(potential_genes)
    obj_div_df = pd.DataFrame() # obj_val: biomass/ growth rate

    if first_n_gene is None: 
        query_gene_list = potential_genes
    elif isinstance(first_n_gene, (int, float)):
        query_gene_list = list(potential_genes)[:first_n_gene]
    else:
        query_gene_list = list(potential_genes)[first_n_gene[0]:first_n_gene[1]]

    for i, current_gene in enumerate(query_gene_list): # iter gene 
#         for i, current_gene in enumerate(obj_flux_df.index[1:][:1]): # iter gene 
        logger.debug('%s in cal', current_gene)
        print(current_gene, 'in cal')
        with model:
            if alpha_df is None: 
                if not search_alpha: search_alpha = 1.02 
                # AlphaFinder class for every model and SG
                AF = MonocultureAlphaFinder(model=model,
                            search_alpha = search_alpha,
                            current_gene = current_gene, 
                            target_obj_val = target_obj_val, 
                            exp_leap=3,
                            precision=precision,
                            acceptance_threshold_upper = acceptance_threshold_upper,
                            acceptance_threshold_lower = acceptance_threshold_lower)
                alpha, obj_value, temp_df = AF.find_feasible_alpha()
            else:
                alpha = alpha_df.loc[current_gene,f'{model.id}_div_opt_alpha'] # get alpha for corresponding gene from alpha_df
                alpha, obj_value,temp_df = Sij_biomass(model, alpha,
                                                    genes = current_gene)
            temp_df['Gene_inhibition'] = current_gene    
            obj_div_df = pd.concat([obj_div_df, temp_df.set_index('Gene_inhibition')],axis=0)
    return obj_div_df

def get_div_obj_df(model_list, target_obj_val, potential_genes=potential_genes, precision=3, detailed_alpha_table=False):
    """Get the objective value for each gene in each model in model_list
    
    Parameters
    ----------
    model_list """
    alpha_obj_df_list = iter_species(model_list, get_single_div_obj_df,
                            target_obj_val=target_obj_val, potential_genes=potential_genes, precision=precision)
    result_df = pd.concat(alpha_obj_df_list, axis=1)  
    if detailed_alpha_table:
        return result_df  
    