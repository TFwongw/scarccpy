import os
import json
import numpy as np

from dataclasses import dataclass, field
from typing import List
import pandas as pd
import concurrent.futures

from scarcc.util import (convert_arg_to_list, rename_columns)
from scarcc.data_analysis.biomass.growth_summary import get_desired_cycle
# from scarcc.sim_engine.simulation_workflow import extract_biomass_flux_df
from .alpha_finder import AlphaFinderConfig

@dataclass(kw_only=True)
class CocultureAlphaFinder(AlphaFinderConfig): # scale normal & knockout to 1-0
    model: List['cobra.Model']
    search_alpha : float
    current_gene : str
    target_obj_val : float = 0.5 
    alpha_table : pd.DataFrame = None
    precision : int = 1 # precision of alpha
#     alpha_table : pd.DataFrame = alpha_table # init AFConfig
    trace_biomass : List = field(default_factory=list)
    trace_obj_val : List = field(default_factory=list) # TODO: combine all trace_XX & convert to dict, follow format of response_record in MonocultureAlphaFinder
    trace_alpha : List = field(default_factory=list)
    gr_ko : float = None
    n_dir : str = ''
    opt_alpha_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    nxt_alpha_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    ko_intercepted : None = None
    acceptance_threshold_lower : float = .95
    acceptance_threshold_upper : float = 1.05
    gr_Normal : float = None
    carbon_source_val : float = 5e-3
    add_nutrient_val : float = 10 # DO NOT set Met limiting in batch culture 
    is_growth_switch = False
    initial_pop : float = 1e-8
    p : None = None
    obj_style: str = 'MAX_OBJECTIVE_MIN_TOTAL'
    data_path = None # replace ./Data/

    # add_nutrient_val : float = field(default_factory=[.08])
    
    def __post_init__(self):
        self.is_monoculture = False
        self.eval_alpha_fun = self.eval_alpha_fun 
        self.get_alpha_use()
        self.E0 = [ele for ele in convert_arg_to_list(self.model) if ele.id == 'E0'][:1]
        # self.gr_Normal = 0.00064644123538125 # ?check effect of lcts conc changed in Div_COMETS
        self.alpha_pair : pd.DataFrame = None
        if not self.exp_leap:
            self.exp_leap = 1.3
        self.E0 = self.model[0]
        self.S0 = self.model[1]
        self.out_fun = self.out_fun
        self.calculate_gr_Normal()
    
    def calculate_gr_Normal(self):
        if not self.gr_Normal:
            if self.current_gene == 'Normal':
                self.gr_Normal = self.calculate_gr_ko()
            else: self.gr_Normal = 0.28672261
        return float(self.gr_Normal)

    def calculate_gr_ko(self):
        print('calculating ko gr')
        self.ko_intercepted = False
        
        if self.current_gene in ['folP', 'folA', 'dadX', 'pgk', 'acnB']:
            self.gr_ko = 0
            return self.gr_ko
        
        # alpha_in_table = self.alpha_table.loc[self.current_gene, self.E0.id] # E0 col
        full_df, *_ =  extract_biomass_flux_df(self.current_gene,n_dir=self.n_dir,alpha_table=self.alpha_table,
                                E0=self.E0, S0=self.S0, mono=False, p=self.p, return_sim=True, ko=True, 
                                carbon_source_val=self.carbon_source_val, add_nutrient_val=self.add_nutrient_val,
                                initial_pop=self.initial_pop, obj_style=self.obj_style)
        
        coculture_biomass_df = full_df.iloc[:,[0]]
        if self.current_gene == 'Normal':
            coculture_biomass_df.to_csv(os.path.join(self.data_path, f'Normal_biomass_{str(self.carbon_source_val)}.csv'))
        
        coculture_biomass_df.columns = rename_columns(coculture_biomass_df)
        
        # desired_cycle = get_desired_cycle(coculture_biomass_df, scale_diff=0.05)
        # # return desired_cycle, full_df, coculture_biomass_df
        # growth_phase = desired_cycle.growth_phase[0]
        # gr_ko = (coculture_biomass_df.loc[growth_phase[1]] - coculture_biomass_df.loc[growth_phase[0]])/(growth_phase[1]-growth_phase[0])
        # gr_ko = float(gr_ko)
        
        gr_ko = self.cal_fitted_gr(coculture_biomass_df)
            
        self.cocdf = coculture_biomass_df 
        # not setting alpha, not upper bound
        self.gr_ko = float(gr_ko)
        
        if self.current_gene == 'Normal':
            print('set_normal')
            self.gr_Normal = gr_ko
        elif self.gr_ko > self.gr_Normal*.9:
            print('intercept')
            self.ko_intercepted = True
            self.opt_alpha_table = self.alpha_table.loc[[self.current_gene]]
        else: print('NNNOT intercept')
        
        print(self.opt_alpha_table)
        
        print(f'self.gr_ko > self.gr_Normal: {self.gr_ko},{self.gr_Normal*.85}')
        print(' gr_ko: ', gr_ko)
        self.trace_biomass.append(full_df.add_suffix(f'_ko'))
        return gr_ko
    
    def get_alpha_use(self):
        if self.current_gene == 'Normal':
            return 0
        Ealpha = self.alpha_table.loc[self.current_gene].iloc[0]
        Salpha = self.alpha_table.loc[self.current_gene].iloc[1]
        self.alpha_use = Ealpha if Ealpha<1e5 else Salpha
        if self.alpha_use > 1e4:
            self.alpha_use = 1.04
            self.exp_leap = 3
        
        if not self.search_alpha and (self.current_gene != 'Normal'):
                self.search_alpha = (self.alpha_use-1)/2+1 if self.alpha_use > 2 else self.alpha_use
        return None
    
    def generate_next_fixed_ratio_table(self):     
        multiplier = (self.search_alpha-1)/(self.alpha_use-1)
        nxt_alpha_table =  (self.alpha_table.loc[[self.current_gene]]-1)*multiplier + 1 # ratio offset by 1, minimum is alpha=1
        # trace alpha table change and corresponding gr
        # self.trace_alpha_table.append(nxt_alpha_table.to_dict())
        self.nxt_alpha_table = nxt_alpha_table
        return nxt_alpha_table

    def cal_fitted_gr(self, coculture_biomass_df):
        def logistic(x, saturation, growth_rate, inflection_point, initial_pop):
            return (saturation / (1 + np.exp(-growth_rate * (x - inflection_point)))) + initial_pop
        
        x = coculture_biomass_df.index
        y = coculture_biomass_df.iloc[:, 0] # E0 column
        # y = coculture_biomass_df # E0 column
        if (float(y.tail(1)) < float(y.head(1))+5e-8):
            growth_rate = 0
            return growth_rate
            
        desired_cycle = get_desired_cycle(coculture_biomass_df, scale_diff=0.05)
        
        growth_phase = desired_cycle.loc[self.current_gene, 'growth_phase']
        bool_growing = desired_cycle.loc[self.current_gene, 'bool_growing']
        c_max_gr = desired_cycle.loc[self.current_gene, 'c_max_gr']
        saturation = float(coculture_biomass_df.loc[growth_phase[1]])
        init_gr = .23
        p0 = [saturation, init_gr, c_max_gr, self.initial_pop]

        print('growing: ',bool_growing)
        if (bool_growing == False) or (bool_growing <1):
            popt, pcov = curve_fit(logistic, x, y, p0=p0)
            print('logi')
            
            if self.current_gene =='Normal':
                print(popt)
            growth_rate = popt[1]
        else:
            print('lolin')
            popt, pcov = np.polyfit(x, np.log(y), 1) # [B,A] -- log(y) = A+Bx
            growth_rate = popt
        if self.current_gene == 'Normal':
            print('popt', popt)
        return growth_rate
        
    def cal_std_gr(self, coculture_biomass_df, gr_Normal):
#         desired_biomass_df = pd.DataFrame(get_cycle_max_gr(coculture_biomass_df), index=['max_gr', 'start', 'end', 'bool_growing'])
        gr = self.cal_fitted_gr(coculture_biomass_df)
        standardized_gr = (float(gr) - self.gr_ko)/(gr_Normal-self.gr_ko) # offset scale self.gr_ko then scale
        if standardized_gr > 1:
            self.exp_leap = 3 # for coculture gr inherently > 1, require large leap
        
        # print('standardized gr: ', standardized_gr)
        self.trace_obj_val.append(standardized_gr)
        print('---GR, stdGR__:',gr, standardized_gr)
        return standardized_gr
    
    def eval_alpha_fun(self, alpha_overwrite=None): # first, search alppha = alpha_table alpha
        if self.ko_intercepted is None:
            _ = self.calculate_gr_ko() # also catch gr_ko ~ gr_Normal- Gene with nonessential reactions
        
        if alpha_overwrite:
            self.search_alpha = alpha_overwrite
        nxt_alpha_table = self.generate_next_fixed_ratio_table()
        # self.trace_alpha_table.append(nxt_alpha_table.to_dict())

        full_df, out_dict, co_sim = extract_biomass_flux_df(self.current_gene,n_dir=self.n_dir,alpha_table=self.nxt_alpha_table,
                                E0=self.E0, S0=self.S0, mono=False, p=p, return_sim=True, ko=False, 
                                carbon_source_val=self.carbon_source_val, add_nutrient_val=self.add_nutrient_val,
                                initial_pop=self.initial_pop, obj_style=self.obj_style)
        # get_BM_df(self.current_gene,n_dir=self.n_dir,alpha_table=nxt_alpha_table,mono=False, return_sim=True)
        self.co_sim = co_sim
        
        self.trace_biomass.append(full_df.add_suffix(f'_{self.i_iter}'))
        
        coculture_biomass_df = full_df.iloc[:,[0]]
        target_gr = self.target_obj_val
        std_gr = self.cal_std_gr(coculture_biomass_df, self.gr_Normal) # ?check this
#         std_gr = self.cal_std_gr(coculture_biomass_df, 1.3)

        self.is_new_ub = (std_gr < target_gr)# sharp bound, although lower estimate, no search for higher alpha if just meet
        self.converge_alpha = self.is_found(self.search_alpha, self.alpha_lb, self.alpha_ub, self.precision)
        
        # self.found_alpha = self.found_alpha or 
        #     (self.std_gr > self.target_obj_val*acceptance_threshold_upper)

        self.obj_found = ((std_gr > self.target_obj_val*self.acceptance_threshold_lower) and 
            (std_gr < self.target_obj_val*self.acceptance_threshold_upper)) 
        
        bool_gr_coculture_exceed_limit = (std_gr > 1)
        self.found_alpha = self.converge_alpha or self.obj_found or self.classify_growth_switch()

        if self.is_new_ub:
            self.nxt_alpha_table = nxt_alpha_table
        print('searc, found, new_un, lu')
        print(self.search_alpha, self.found_alpha, self.is_new_ub, self.alpha_lb, self.alpha_ub)
        print('convg, obj_f, exceed, gr_sw')
        print(self.converge_alpha, self.obj_found, bool_gr_coculture_exceed_limit, self.classify_growth_switch())
        # update optimal df
        
        obj_req = ((std_gr > self.target_obj_val*0.95) and (std_gr < self.target_obj_val*1.1)  # harsh lb
                    or self.is_growth_switch)
        
        temp_df = pd.DataFrame.from_dict(
            {self.current_gene :
                {'std_growth_rate' : std_gr,
                'lb_ub' : (self.alpha_lb, self.alpha_ub),
                'converge_alpha' : self.converge_alpha,
                'obj_found': self.obj_found,
                'bool_gr_coculture_exceed_limit': bool_gr_coculture_exceed_limit,
                'classify_growth_switch': self.is_growth_switch,
                'Non_essential': self.ko_intercepted
                }}, 
                orient='index')
        nxt_alpha_table = pd.concat([nxt_alpha_table, temp_df], axis=1)
        nxt_alpha_table.index.name='Gene_inhibition'
        self.trace_alpha.append(nxt_alpha_table.to_dict())
        
        
        if obj_req or (self.opt_df is None) or self.converge_alpha: # store only qualified alpha
            if self.opt_df is not None:
                nxt_alpha_table['obj_req_satisfied'] = True
            self.opt_df = [full_df, nxt_alpha_table]
            # self.opt_alpha_table = pd.concat([nxt_alpha_table, temp_df], axis=1)
            self.opt_alpha_table = nxt_alpha_table
        return full_df 
    
    def out_fun(self):
        result = {}
        for i in range(len(self.trace_obj_val)):
            result[f'iter_{i+1}'] = {'obj_val': self.trace_obj_val[i], 'alpha_table': self.trace_alpha[i]}
        
        biomass_df = pd.concat(self.trace_biomass, axis=1) if self.trace_biomass else pd.DataFrame()
        return biomass_df, self.opt_alpha_table, {self.current_gene : result}
#         return 'Temp Done' if (n_iter<=10 or n_iter ==2) else 'end at: '+ str(n_iter)

##--------Perform Coculture Secrch---------##

def coculture_search_job(**kwargs):    
    AF = CocultureAlphaFinder(**kwargs)
    AF.calculate_gr_ko()
    return AF.find_feasible_alpha()

def run_coculture_search_mp(potential_genes, data_path, filename, n_processor, **kwargs):
    def save_trace_biomass(trace_biomass):
        trace_biomass = pd.concat(trace_biomass,axis=1) # alpha_table
        trace_biomass.columns = rename_columns(trace_biomass)
        trace_biomass.to_csv(os.path.join(data_path, f'biomass_{filename}.csv'))
    
    result_list = list()
    with concurrent.futures.ProcessPoolExecutor(n_processor) as executor:
        for i, current_gene in enumerate(potential_genes):
            future = executor.submit(coculture_search_job, current_gene=current_gene, n_dir=i, **kwargs)
            result_list.append(future)
    
    trace_biomass, opt_alpha_list, result_dict_list = zip(*[r.result() for r in result_list])
    opt_alpha = pd.concat(opt_alpha_list) # alpha_table
        
    opt_alpha.columns = rename_columns(opt_alpha)
    opt_alpha.to_csv(os.path.join(data_path, f'alpha_table_{filename}.csv'))
    print(opt_alpha, result_dict_list)
    
    trace_dict = {k: v for d in result_dict_list for k, v in d.items()}
    with open(os.path.join(data_path, f"trace_record_{filename}.json", "w")) as outfile:
        json.dump(trace_dict, outfile)

    save_trace_biomass(trace_biomass)
    print('finished alpha search')
    return opt_alpha_list, trace_dict