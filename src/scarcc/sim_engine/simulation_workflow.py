import os
import logging
import pandas as pd
from typing import List
from dataclasses import dataclass, field
from collections import defaultdict
import itertools

import concurrent.futures
import cometspy as c

from scarcc.utils import convert_arg_to_list
from scarcc.data_analysis.growth.growth_rate import MethodDataFiller

from .simulation_configuration import LayoutConfig
from scarcc.preparation.perturbation import get_alphas_from_tab, alter_Sij
from .flux_extraction import extract_biomass_flux_df

logger = logging.getLogger(__name__)

def sim_culture(layout, p=None, base = None):
    # flux data needs biomass to determine the ~, only return sim object for the outpout object
    # separate function into mono & coculture to prevent using wrong layer
    if isinstance(layout, list):
        if len(layout) > 1:
            raise ValueError("The list 'layout' should contain only one element.")
        (layout,) = layout # one element unpacking from iter_species
    
    sim = c.comets(layout, p)
    sim.working_dir = os.path.join(base, '') # make sure it is a directory instead of file
    print(sim.working_dir)
    
    try:
        sim.run()
    except:
        logging.exception(f"{sim.run_output}")
        print(f"{sim.run_output}")
    return sim

@dataclass(kw_only=True)
class SimulateCombinedAntibiotics(LayoutConfig):
    current_gene: str
    p: 'comets.p'
    alpha_table: str
    base: str = None # base as __file__ or working directory
    checker_suffix: str = None
    return_sim: bool = False
    ko: bool = False

    # default values for output
    working_dir: str = None # passed to comets object
    biomass_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    co_sim_object: 'comets.simulation' = field(default_factory=list)
    mono_sim_object_list: List['comets.simulation'] = field(default_factory=list)
    sim_object_list: List['comets.simulation'] = field(default_factory=list)

    def __post_init__(self):
        super().__post_init__()
        self.current_gene = convert_arg_to_list(self.current_gene)

        path_elements = [self.base, 'SimChamber', '.'.join(self.current_gene)] if 'SimChamber' not in self.base else [self.base, '.'.join(self.current_gene)]
        self.working_dir = os.path.join(*path_elements) # '' as specification of directory where COMETS files are stored
        os.makedirs(self.working_dir, exist_ok=True)
        
        # filepath
    def cleanup(self):
        try:
            os.rmdir(self.working_dir)
        except:
            print(f'Failed to remove {self.working_dir}, needs remove files manually')

    def get_BM_df(self):
        with self.E0 as m_E0, self.S0 as m_S0:
            metabolic_model_list = [m_E0, m_S0]
            if 'Normal' not in self.current_gene:
                alphas = [get_alphas_from_tab(model, genes=self.current_gene, alpha_table=self.alpha_table) for model in metabolic_model_list]
                _ = [alter_Sij(model, alphas=alpha, genes=self.current_gene, ko=self.ko) for model, alpha in zip(metabolic_model_list, alphas)]

            E_model, S_model = self.set_comets_model()
            co_layout, E0_layout, S0_layout = self.set_layout_object()

            if self.co:
                logger.debug(f'{self.p.all_params["maxCycles"]} co_p')
                self.co_sim_object = sim_culture(self.co_layout, p=self.p, base=self.working_dir)
                self.sim_object_list.append(self.co_sim_object)

            if self.mono:
                monoculture_to_run = dict()
                if self.mono_E:
                    monoculture_to_run[self.E_model] = self.monoE_layout
                if self.mono_S:
                    monoculture_to_run[self.S_model] = self.monoS_layout
                self.mono_sim_object_list = [sim_culture(layout, p=self.p, base=self.working_dir) for layout in monoculture_to_run.values()]
                self.sim_object_list.extend(self.mono_sim_object_list)
            
        self.biomass_df, self.flux_df = extract_biomass_flux_df(self.E0, self.S0, self.sim_object_list, alpha_table=self.alpha_table, current_gene=self.current_gene)
        self.cleanup()
        return self.biomass_df, self.flux_df

def read_alpha_table(data_directory, alpha_table_suffix):
    # check file exist
    # ? check alpht_table contains all SG
    file_dir = os.path.join(data_directory, f'alpha_table_{alpha_table_suffix}.csv')
    if os.path.isfile(file_dir):
        alpha_table = pd.read_csv(os.path.join(data_directory, f'alpha_table_{alpha_table_suffix}.csv'))
        return alpha_table
    else:
        # ? run procedure for creating alpha_table
        raise FileNotFoundError(f'File {file_dir} does not exist')

def nested_dict(d):
    result = {}
    for keys, value in d.items():
        current_dict = result
        for key in keys[:-1]:
            current_dict = current_dict.setdefault(key, {})
        current_dict[keys[-1]] = value
    return result

def unpack_future_result_per_key(result_list):
    keys = ['biomass', 'flux']
    # return {k: v for k, v in zip(keys, zip(*[r.result() for r in result_list]))}
    return dict(zip(keys, zip(*result_list)))

def concat_result(result_dict):
    def concat_df(key, df_list):
        if 'biomass' in key:
            return pd.concat(df_list, axis=1)
        return pd.concat(df_list, axis=0)
    return {k: concat_df(k, v) for k, v in result_dict.items()}

