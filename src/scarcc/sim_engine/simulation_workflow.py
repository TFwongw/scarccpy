import os
import logging
import pandas as pd
from typing import List
from dataclasses import dataclass, field
# from collections import defaultdict
import itertools
import concurrent.futures
import cometspy as c

from scarcc.utils import convert_arg_to_list
from scarcc.preparation.perturbation import get_alphas_from_tab, alter_Sij
from scarcc.sim_engine.output import MethodDataFiller
from .simulation_configuration import LayoutConfig
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
        print(path_elements)
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
        alpha_table = pd.read_csv(os.path.join(data_directory, f'alpha_table_{alpha_table_suffix}.csv'), index_col=0)
        return alpha_table
    else:
        # ? run procedure for creating alpha_table
        raise FileNotFoundError(f'File {file_dir} does not exist')

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

def check_files_exist(data_directory, file_paths):
    missing_files = []
    short_file_paths = file_paths.copy()
    file_paths = [os.path.join(data_directory, file_path) for file_path in file_paths]
    for short_file_path, file_path in zip(short_file_paths, file_paths):
        if not os.path.exists(file_path):
            missing_files.append(short_file_path)
    if missing_files:
        if any('SG' in file_path for file_path in missing_files):
            print('Required files for SG are missing for calculation of DG data frames, please check the file paths or set generate_SG_list to True in run_sim_workflow')
        raise FileNotFoundError(f'Files {missing_files} do not exist in {data_directory}')

def check_SG_DG_format(SG_list, DG_list, generate_SG_list):
    DG_format_list = True
    if DG_list is not None:
        if len(DG_list) == 2:
            if isinstance(DG_list[0], str):
                DG_format_list = False
                if '.' not in DG_list[0]:
                    raise ValueError('DG_list should be a list of lists, each list contains gene names')
    if generate_SG_list:
        if not DG_format_list:
            DG_list = [ele.split('.') for ele in DG_list]
        new_SG_list = itertools.chain(*DG_list)
        SG_list = list(set(SG_list).union(set(new_SG_list)))
    if 'Normal' not in SG_list:
        SG_list = ['Normal'] + SG_list
    return SG_list, DG_list

def get_simulation_output(**kwargs):
    simulation = SimulateCombinedAntibiotics(**kwargs)
    return simulation.get_BM_df()

def run_sim_workflow(method_list, data_directory, SG_list=None, DG_list=None, generate_SG_list=False, max_cpus=12, **kwargs):
    SG_list, DG_list = check_SG_DG_format(SG_list, DG_list, generate_SG_list) 
    available_cpus = min(os.cpu_count(), max_cpus)
    method_list = [ele.replace('alpha_table_', '') for ele in method_list]

    required_files = [f'alpha_table_{method}.csv' for method in method_list]
    gene_list_dict = {}
    for XG_list, XG in zip([SG_list, DG_list], ['SG', 'DG']):
        if XG_list is not None:
            gene_list_dict[XG] = XG_list
        else:
            if XG == 'SG':
                required_files.append(f'BM_{XG}_{method}.csv' for method in method_list)
    check_files_exist(data_directory, required_files)
    alpha_table_dict = {method: read_alpha_table(data_directory, method) for method in method_list}

    task_dict = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=available_cpus) as executor:
        for method, XG in itertools.product(method_list, gene_list_dict.keys()):
            # replace with function for running sim
            task_dict[method, XG] = [executor.submit(get_simulation_output, current_gene=current_gene, alpha_table=alpha_table_dict[method], **kwargs)
                                                for current_gene in gene_list_dict[XG]]
    df_container = {key: [future.result() for future in future_list]
                    for key, future_list in task_dict.items()}
    df_container = {key: unpack_future_result_per_key(result_list) for key, result_list in df_container.items()} # list of [biomass, flux] into {biomass : biomass_list, flux : flux_list}
    df_container = {k: concat_result(sub_container) for k, sub_container in df_container.items()} # column-wise for biomass, row-wise for flux data frame concatenation

    mdf = MethodDataFiller(df_container, data_directory)
    mdf.fill_container()
    mdf.write_to_csv()
    return df_container