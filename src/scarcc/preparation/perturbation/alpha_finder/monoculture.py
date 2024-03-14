"""Search for IC values by percentage of biomass production in FBA under perturbation"""

from dataclasses import dataclass, field
import logging
import pandas as pd
from typing import List, Dict

from scarcc.preparation.perturbation import alter_Sij, iter_species
from scarcc.utils import convert_arg_to_list
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
    """FBA monoculture alpha finder to restrict biomass production at Normalized target biomass
    *args and **kwargs are passed to the AlphaFinderConfig parent class

    Parameters
    ----------
    model : 'cobra.Model'
        Single cobra model
    exp_leap : float
        Exponential search leap
    is_monoculture : bool
        Whether the alpha search is done in FBA, monoculture
    trace_obj_val : List
        List of objective values result from the search process
    norm_obj_val : float
        Normalized objective value
    response_record : Dict
        Record flux and other summary information under different alpha from the search process
    opt_df : pd.DataFrame
        Detailed data frame to be returned with flux and reaction information for the optimal alpha  

    Attributes
    ----------
    response_record : Dictionary with the following structure
        {current_gene: {model: {'Normalized target biomass': target_obj_val, 
                                'response': {
                                    obj_val_interval: {
                                        'search_alpha' : search_alpha,
                                        'precise_obj_val' : precise_obj_val, 
                                        'interval_abs_diff' : interval_abs_diff,
                                        'target_abs_diff' : target_abs_diff
                                    }
                                }}}}
    
    """
    # TODO: Implementation using linear programming to skip search process

    model : 'cobra.Model' # single cobra model
    exp_leap : float = 2
    is_monoculture : bool = True
    trace_obj_val : List = field(default_factory=list)
    norm_obj_val : float = None # gr_Nomral in coculture case
    response_record: Dict = None # for iterative update in checkerboard
    opt_df : pd.DataFrame = None

    def __post_init__(self):
        if not self.acceptance_threshold_lower:
            self.acceptance_threshold_lower = .85
        if not self.acceptance_threshold_upper:
            self.acceptance_threshold_upper = 1.1
        self.norm_obj_val = self.model.slim_optimize()
        if not self.response_record:
            self.response_record = {self.current_gene: {self.model: {'Normalized target biomass': convert_arg_to_list(self.target_obj_val), 'response': {}}}}

    def fill_response_record(self, obj_val, is_lowest_abs_diff=False):
        """Fill response record with the objective value, if objective value falls into previously visited category, update the record

        Parameters
        ----------
        obj_val : float
            Objective value
        is_lowest_abs_diff : bool
            Whether the objective value has the lowest absolute difference from the target objective value
        """
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
        """Evaluate the current alpha and update the optimal dataframe"""
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
        """Output the result of the alpha search
        
        Returns
        -------
        Tuple
            Series of alpha, biomass, and the optimal dataframe"""
        opt_df = self.opt_df
        opt_df['is_growth_switch'] = self.classify_growth_switch()
        logger.debug('Gene: %s with alpha %s and obj_val %s', self.current_gene, opt_df['alpha'][0], opt_df[f'biomass'][0])
        opt_df.insert(loc=2, column='Normalized_biomass', value=opt_df['biomass']/self.norm_obj_val)
        opt_df.columns = [f'{self.model.id}_'+element for element in opt_df.columns]
        self.opt_df = opt_df
        return (opt_df[f'{self.model.id}_alpha'], opt_df[f'{self.model.id}_biomass'], opt_df) # alpha_feasible, and upper bound of alpha

def get_summary_df(model, alpha, obj_val, rct_ids = 'DHFR', sol_sol = None) -> pd.DataFrame: 
    """Dataframe summary from the optimization result

    Parameters
    ----------
    model : cobra.Model
    alpha : float
    obj_val : float
        target biomass
    rct_ids : str
        reaction ids
    sol_sol : FBA solution, Optional - flux distribution summary

    Returns
    -------
    pd.DataFrame with the following columns
        - 'alpha' (float): optimal searched alpha
        - 'biomass' (float): biomass yield associated with the perturbation characterized by alpha  
        - 'forward_reactions_id' (str): reaction ids in forward direction
        - 'forward_reactions' (str): reaction names in forward direction
        - 'forward_flux_values' (list): flux values in forward direction
        - 'opposite_reactions' (str): reaction ids in opposite direction
        - 'opposite_flux_values' (list): flux values in opposite direction
        - 'Net_Flux_Boolean' (str): 'Zero Flux' if net flux is zero, 'Net Flux' otherwise
        - 'Net_Flux' (float): net flux
    """
    if type(alpha) != int and type(alpha) != float:
        alpha = str(alpha)

    if sol_sol is None:
        return pd.DataFrame([{'alpha': alpha,'biomass': obj_val}])
    
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
    net_flux = sum(abs(element) for element in set(flux_dir) | set(flux_rev))
    net_flux_bool = 'Zero Flux' if net_flux ==0 else 'Net Flux'

    summary_df = pd.DataFrame({f'alpha': alpha,
                                f'biomass': obj_val,
                                f'forward_reactions_id': ', '.join(rct_ids),
                                f'forward_reactions': ' '.join(rct_dir),                                    
                                f'forward_flux_values':[flux_dir],
                                f'opposite_reactions': ', '.join(rct_rev),     
                                f'opposite_flux_values':[flux_rev],
                                f'Net_Flux_Boolean':[net_flux_bool],
                                f'Net_Flux':[net_flux]}
                                )
    return(summary_df)

def Sij_biomass(model: "cobra.Model", alphas: float = 1, genes: str = 'folA', slim_opt: bool = False):
    """Get biomass form FBA optimization with given gene and alpha
    
    Parameters
    ----------
    model : cobra.Model
    alphas : float
    genes : str
    slim_opt : bool
        If True, use slim_optimize
        If False, use optimize
        
    Returns
    -------
    float if slim_opt is True
        Objective value
    list(alpha, biomass, summary_df) if slim_opt is False
    """
    with model:
        rct_ids = alter_Sij(model, alphas, genes) # same

        if (slim_opt == False):
            sol_sol = model.optimize()
            obj_val = sol_sol.objective_value 
            summary_df = get_summary_df(model, alphas, obj_val, rct_ids, sol_sol.fluxes)  
            return(alphas, obj_val, summary_df) 
        else:
#             return(pd.DataFrame([{f'biomass': model.slim_optimize()}]))
            return(model.slim_optimize())

def get_single_div_biomass_df(model: "cobra.Model", target_obj_val: float, potential_genes: list, precision: int = 3, search_alpha: float= 1.02,
                        alpha_table: pd.DataFrame = None, acceptance_threshold_upper: float = 1.1, acceptance_threshold_lower: float=.95) -> pd.DataFrame: 
    """Get the objective value for each gene in a single model
    If alpha_table is not provided, perform a search for the optimal alpha 
    and record the corresponding objective value
    
    Parameters
    ----------
    model : cobra.Model
    target_obj_val : float
        The desired target objective value
    potential_genes : list
        List of potential genes to generate alpha_table
    precision : int
        The precision of the alpha
    search_alpha : float
        initial alpha to start search
    alpha_table : pd.DataFrame, Optional
        The alpha table
    acceptance_threshold_upper : float
        The upper acceptance threshold for the target objective value
    acceptance_threshold_lower : float
        The lower acceptance threshold for the target objective value

    Returns
    -------
    pd.DataFrame
        summary_df        
    """
    alpha_biomass_df = pd.DataFrame() # obj_val: biomass/ growth rate

    for current_gene in potential_genes: # iter gene 
        logger.debug('Searching Model: %s, Gene: %s ', model.id, current_gene)
        if alpha_table is None: 
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
            alpha = alpha_table.loc[current_gene,f'{model.id}_alpha'] # get alpha for corresponding gene from alpha_df
            alpha, obj_value, temp_df = Sij_biomass(model, alpha,
                                                genes = current_gene)
        temp_df['Gene_inhibition'] = current_gene 
        alpha_biomass_df = pd.concat([alpha_biomass_df, temp_df.set_index('Gene_inhibition')],axis=0)
    return alpha_biomass_df

def get_alpha_biomass_df(model_list: List["cobra.Model"]=None, target_normalized_biomass: float=None, potential_genes: List=potential_genes,
                    precision: float=3, detailed_alpha_table: bool=False, alpha_biomass_df: pd.DataFrame=None):
    """Get the objective value for each gene in each model in model_list
    
    Parameters
    ----------
    model_list : List
        List of cobra models
    target_normalized_biomass : float
        The desired target objective value
    potential_genes : List
        List of potential genes 
    precision : float
        The precision of the alpha
    detailed_alpha_table : bool
        boolean to return detailed alpha table or not
    alpha_biomass_df : pd.DataFrame, Optional
    
    Returns
    -------
    pd.DataFrame
        summary_df if detailed_alpha_table is False
        alpha and normalized_biomass columns of summary_df otherwise
    """
    if alpha_biomass_df is None:
        alpha_obj_df_list = iter_species(model_list, get_single_div_biomass_df,
                                target_obj_val=target_normalized_biomass, potential_genes=potential_genes, precision=precision)
        result_df = pd.concat(alpha_obj_df_list, axis=1)
    else:
        result_df = alpha_biomass_df
    if detailed_alpha_table:
        return result_df
    return result_df.filter(regex='alpha|Normalized_biomass')
    