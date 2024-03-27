"""Module for alpha finder base class"""

from dataclasses import dataclass
from typing import List
from abc import ABC, abstractmethod
from typing import Union
import logging

import cobra

logger = logging.getLogger(__name__)

@dataclass(kw_only=True)
class AlphaFinderConfig(ABC):
    """Configure the alpha finder class with procedures update next alpha and classification of gene

    Parameters
    ----------
    current_gene : str
        gene name
    trace_obj_val : dict
        objective value of the gene
    is_growth_switch : bool
        whether the gene is growth switch, biomass drop abruptly within certain precision of alpha, mostly transport reaction in energy pathway
    early_stop : bool
        reached maximum iteration - without reaching the target precision of alpha
    ko_intercepted : bool
        stop finding alpha since the gene is non-essential, producing no effect to change in biomass
    is_monoculture : bool
        whether the alpha search is done in FBA, monoculture
    model : Union["cobra.Model", List["cobra.Model"]]
        single species in FBA
        list of species in COMETS coculture simulation 
    ever_eval : bool
        initialize first run in coculture
    target_obj_val : float
        target biomass value produced by the perturbation by alpha
    search_alpha : float
        alpha value to be started with, updated in each iteration
    alpha_lb : float
        lower bound of alpha, objective value is higher than target one, updated in each iteration
    alpha_ub : float
        upper bound of alpha, objective value is lower than target one, updated in each iteration
    precision : int
        precision of alpha
    is_new_ub : bool
        whether the upper bound is updated in the current iteration
    exp_leap : float
        multiplier factor of search_alpha in exponential search if new_ub is False, implement binary search otherwise
    found_alpha : bool
        whether the alpha in inhibiting biomass yield to target biomass is found
    iter_max : int
        maximum iteration for alpha search
    i_iter : int
        current iteration of alpha search, initialized with 0 
    acceptance_threshold_upper : float
        multiplier of acceptance threshold for target biomass, value higher than 1
    acceptance_threshold_lower : float
        lower bound of acceptance threshold for target biomass, value between 0 and 1

    Remarks
    -------
    Configuration of alpha finder class for MonocultureAlphaFinder and CocultureAlphaFinder, DO NOT use this class directly
    """

    # gene characteristic and records
    current_gene: str = None
    trace_obj_val: dict = None
    is_growth_switch: bool = None
    early_stop: bool = False
    ko_intercepted: bool = False

    # culture identification
    is_monoculture: bool = None
    model: Union["cobra.Model", List["cobra.Model"]] = None # single model for monoculture, list for coculture
    ever_eval: bool =  False # initialize first run coculture

    # alpha search related
    target_obj_val : float = None
    search_alpha: float = None
    alpha_lb : float = 1+1e-7
    alpha_ub : float = 1e5
    precision : int = 2 # precision of alpha
    is_new_ub : bool = None
    exp_leap : float = 2
    found_alpha : bool = None
    iter_max: int = 25
    i_iter: int = 0

    # evaluation alpha choice
    acceptance_threshold_upper : float = None
    acceptance_threshold_lower : float = None

    @abstractmethod
    def eval_alpha_fun(self):
        pass

    @abstractmethod
    def out_fun(self):
        pass

    @staticmethod
    def get_next_alpha(search_alpha, alpha_lb, alpha_ub, is_new_ub, exp_leap=2):
        """Update alpha search value in each iteration if target biomass is not reached"""
        if is_new_ub is False:  # raise lb, search higher alpha
            alpha_lb = search_alpha
            if (search_alpha*exp_leap <alpha_ub and is_new_ub is False):
                search_alpha *= exp_leap
            else: # search alpha not get updated, then binary search
                search_alpha = (search_alpha+alpha_ub)/2
        else: # search alpha not accepted, lower ub, lower search_alpha
            alpha_ub = search_alpha
            search_alpha = max((search_alpha+alpha_lb)/2, 1+1e-5) # terminate at 1+1e-5
        return search_alpha, alpha_lb, alpha_ub
    
    @staticmethod
    def is_found(search_alpha, alpha_lb, alpha_ub, precision):
        """Check whether the alpha is found within the precision"""
        if (search_alpha<2) and (precision <3):
            precision = 3
        if alpha_lb > 15:
            precision = 1
        return (round(search_alpha,precision)==round(alpha_ub,precision) or
                round(search_alpha,precision)==round(alpha_lb,precision) or
                (search_alpha>9e4))

    def classify_growth_switch(self, min_attain=.9, max_attain=.3) -> bool:
        """Classify the gene as growth switch
        biomass drop abruptly within certain precision of alpha, mostly transport reaction in energy pathway
        need to reach min_attain and max_attain for biomass value to classify as growth switch, 
        risk of misclassification otherwise
        Precision adjusted for alpha range

        Parameters
        ----------
        min_attain : float
            minimum biomass value to be attained
        max_attain : float
            maximum biomass value to be attained

        Returns
        -------
        bool
            whether the gene is growth switch
        """
        alpha_range = (self.alpha_ub - self.alpha_lb)
        alpha_range_narrow = (((alpha_range < 0.15) and (self.alpha_ub < 3)) or # circumvent precision <2 in subsequent evaluation
                                ((alpha_range < .3) and (5<self.alpha_lb<10)) or
                                ((alpha_range < 1.3) and (self.alpha_lb > 10)))
        alpha_range_req =  (alpha_range_narrow
                            and (self.alpha_lb > 1.001)) # ever update lb and ub
        
        obj_req = (not any(0.3 < val < 0.8 for val in self.trace_obj_val) and # mrdA forever > 1 for low dose
                    (min(self.trace_obj_val) < min_attain) and
                    (max(self.trace_obj_val) > max_attain))
        
        self.is_growth_switch = (alpha_range_req and obj_req) # or (self.alpha_ub < 1.018) purT
        if self.net_flux is not None and self.net_flux < 1e-8: # already without flux -> non-essential
                self.is_growth_switch = False
        return self.is_growth_switch
    
    def find_feasible_alpha(self):
        """If non-essential, producing no effect to change in biomass, jump to output function
        If essential, find the alpha that can inhibit the biomass yield to target biomass

        Returns
        -------
        output of out_fun defined in child class
        """
        def eval_continue():
            # see if conditions for stopping are met
            if self.is_monoculture: # as culture_flag
                return not self.found_alpha
            if self.i_iter>self.iter_max:
                self.early_stop = True
                logger.debug('Early stopped coculture')
            return ((not self.found_alpha) or (self.i_iter<2)
                    or (self.alpha_lb < 1.01 and 1.018< self.alpha_ub)) # force iter -ensure optimal


        if self.ko_intercepted: # req run ko_gr first, otherwise cannot catch
            logger.debug(f'Intercept Non-essential: {self.current_gene}, calculating with current alpha')
            return self.out_fun() # sim_culture with current alpha_table already ran
                
        if not self.ever_eval:
            self.eval_alpha_fun() # initial ko is without alpha, only used for identify gr_ko

        # stop = not (self.alpha_lb < self.search_alpha < self.alpha_ub) and 
        
        # while eval_continue() and not stop: 
        while eval_continue(): 
            self.search_alpha, self.alpha_lb, self.alpha_ub = self.get_next_alpha(
                self.search_alpha, self.alpha_lb, self.alpha_ub, self.is_new_ub,
                self.exp_leap)
            self.eval_alpha_fun()
            self.i_iter+=1
        if (not self.is_monoculture):
            logger.debug(f'Stopped at iter {self.i_iter}') if (self.i_iter<=self.iter_max or self.i_iter ==2) else logger.debug('Success search, end at: ', str(self.i_iter))
        
        return self.out_fun()
