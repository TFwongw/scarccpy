import pandas as pd
import os
import cometspy as c
import logging
import concurrent.futures
import itertools

from scarcc.preparation.find_directory import find_directory
from scarcc.preparation.metabolic_model import BasicModel

from scarcc.sim_engine.output import MethodDataFiller
from scarcc.sim_engine.simulation_workflow import (SimulateCombinedAntibiotics,
                                                    read_alpha_table,
                                                    unpack_future_result_per_key,
                                                    concat_result, run_sim_workflow)

logging.basicConfig(filename='test_simulation_workflow.log',
                    level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')

method_list = ['m1']
SG_list = ['folA', 'folP']
DG_list = [['folA', 'folP']]

data_directory = find_directory('Data', os.path.abspath(''))
model_directory = find_directory('models', os.path.abspath(''))
sim_chamber_directory = find_directory('SimChamber', os.path.abspath(''))

E0, S0, all_components = BasicModel(flux_weighting=True, model_directory=model_directory).load_ES_models()

p = c.params()
p.set_param("defaultKm", 0.00001) # M
p.set_param("defaultVmax", 10) #mmol/gDw/hr
p.set_param("maxCycles", 30)
# p.set_param("maxCycles", 150)
p.set_param("timeStep", 1)
p.set_param('writeFluxLog', True)


simulation_kwargs = {
    'E0': E0,
    'S0': S0,
    'base': sim_chamber_directory,
    'p': p}

if __name__ == '__main__':
    df_container = run_sim_workflow(method_list=['m1'], data_directory=data_directory,
                                    base = sim_chamber_directory, SG_list=['folA', 'folP'], DG_list=[['folA', 'folP']], **simulation_kwargs)
    
    assert ('m1', 'SG') in df_container
    assert 'biomass' in df_container[('m1', 'SG')]
    assert 'flux' in df_container[('m1', 'SG')]
    assert 'normalized_growth_rate' in df_container[('m1', 'SG')]
    assert 'drug_response_classification' not in df_container[('m1', 'SG')]
    assert 'drug_response_classification' in df_container[('m1', 'DG')]