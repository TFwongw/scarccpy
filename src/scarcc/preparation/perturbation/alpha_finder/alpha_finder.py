from dataclasses import dataclass
import pandas as pd
from typing import List

@dataclass(kw_only=True)
class AlphaFinderConfig:
    alpha_table = None
    opt_df : pd.DataFrame = None # final output
    alpha_lb : float = 1+1e-7
    alpha_ub : float = 1e5 
    found_alpha : bool = None
    is_new_ub : bool = None
    exp_leap : float = 2
    target_obj_val : dict = None  
    ko_intercepted : bool = None
    ever_eval =  False
    iter_max = 25
    i_iter = 0    

    def __post_init__(self): # not post_inited
        self.is_monoculture = None
        self.out_fun = None
        self.eval_alpha_fun = None
    
    @staticmethod
    def get_next_alpha(search_alpha, alpha_lb, alpha_ub, is_new_ub, exp_leap=2, is_monoculture=True):
        # print('gna', search_alpha, alpha_lb, alpha_ub, is_new_ub)
        if (is_new_ub == False):  # raise lb, search higher alpha
            alpha_lb = search_alpha
            # if ((search_alpha*exp_leap <alpha_ub) and is_monoculture): # exponential search(large step)
            # if ((alpha_ub > 1e4) and (alpha_lb> 1e3)): # when alpha ub not get updated, exponential search(large step)
            if (search_alpha*exp_leap <alpha_ub and is_new_ub == False):
                search_alpha *= exp_leap
            else: # search alpha not get updated, then binary search
                search_alpha = (search_alpha+alpha_ub)/2  
        else: # search alpha not accepted, lower ub, lower search_alpha
            # start binary search, in (alpha_lb, search_alpha)
            alpha_ub = search_alpha
            search_alpha = max((search_alpha+alpha_lb)/2, 1+1e-5) # terminate at 1+1e-5
        # print(search_alpha, alpha_lb, alpha_ub)
        return search_alpha, alpha_lb, alpha_ub
    
    @staticmethod
    def is_found(search_alpha, alpha_lb, alpha_ub, precision):
        if (search_alpha<2) and (precision <3):
            precision = 3
        if alpha_lb > 15:
            precision = 1
        return (round(search_alpha,precision)==round(alpha_ub,precision) or
                round(search_alpha,precision)==round(alpha_lb,precision) or
                (search_alpha>9e4))

    def classify_growth_switch(self, min_attain=.9, max_attain=.3):
        # eval all true
        alpha_range = (self.alpha_ub - self.alpha_lb)
        alpha_range_narrow = (((alpha_range < 0.15) and (self.alpha_ub < 3)) or # circumvent precision <2 in subsequent evaluation
                                ((alpha_range < .3) and (5<self.alpha_lb<10)) or
                                ((alpha_range < 1.3) and (self.alpha_lb > 10)))
        alpha_range_req =  (alpha_range_narrow # 
                            and (self.alpha_lb > 1.01)) # ever update lb and ub
        
        obj_req = (not any(0.3 < val < 0.8 for val in self.trace_obj_val) and # mrdA forever > 1 for low dose
                    (min(self.trace_obj_val) < min_attain) and 
                    (max(self.trace_obj_val) > max_attain))
        
        if alpha_range_narrow and all(val >1 for val in self.trace_obj_val): # small dose all > Normal
            obj_req = True
                
        # ? req more evalluation for alpha inside the bound?
        self.is_growth_switch = (alpha_range_req and obj_req) or (self.alpha_ub < 1.018)
        return self.is_growth_switch
    
    def find_feasible_alpha(self): 
        print(self.current_gene, self.ko_intercepted)
        def eval_continue():
            if self.is_monoculture: # as culture_flag
                return not self.found_alpha
            if self.i_iter>self.iter_max:
                self.early_stop = True
                print('Early stopped coculture')
            return ((not self.found_alpha) or (self.i_iter<2) 
                    or (self.alpha_lb < 1.01 and 1.018< self.alpha_ub)) # force iter -ensure optimal

        if self.ko_intercepted: # req run ko_gr first, otherwise cannot catch
            print(f'Intercept Non-essential: {self.current_gene}, calculating with current alpha')
            return self.out_fun() # sim_culture with current alpha_table already ran
                
        if not self.ever_eval:
            self.eval_alpha_fun() # initial ko is without alpha, only used for identify gr_ko

        stop = not (self.alpha_lb < self.search_alpha < self.alpha_ub)
        
        while eval_continue() and not stop: 
            self.search_alpha, self.alpha_lb, self.alpha_ub = self.get_next_alpha(
                self.search_alpha, self.alpha_lb, self.alpha_ub,  self.is_new_ub, 
                self.exp_leap, is_monoculture=self.is_monoculture)
            self.eval_alpha_fun()
            self.i_iter+=1
        if (not self.is_monoculture):
            print(f'Stopped at iter {self.i_iter}') if (self.i_iter<=self.iter_max or self.i_iter ==2) else print('Success search, end at: ', str(self.i_iter))
        
        return self.out_fun()

potential_genes = [
    'glyA', 'gltA', 'tktA', 'dapF', 'dapB', 'acnB', 'pgk', 'talB', 'mrcA', 'pyrE', 'dapD', 'pfkA', 'gdhA', 'folA', 'mrdA', 'thrB',
    'dapA', 'serC', 'argD', 'thrC', 'aceF', 'pykF', 'dadX', 'folC', 'pyrD', 'trpA', 'serB', 'fbp', 'eno', 'pgi', 'pheA',
    'gcvH', 'gnd', 'murA', 'aroA', 'guaB', 'glnA', 'yrbG', 'folP', 'purU', 'serA', 'gltD', 'purT', 'ackA', 'purN', 'rffG',
    'gapA'
]