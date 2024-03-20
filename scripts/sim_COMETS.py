import pandas as pd
import os
import cometspy as c
import logging

from scarcc.preparation.find_directory import find_directory
from scarcc.preparation.metabolic_model import BasicModel

def main():
    logging.basicConfig(filename='scarcc_simulation.log',
                        level=logging.DEBUG,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filemode='w')

    data_directory = find_directory('Data', os.path.abspath(''))
    model_directory = find_directory('models', os.path.abspath(''))
    sim_chamber_directory = find_directory('SimChamber', os.path.abspath(''))

    E0, S0, all_components = BasicModel(flux_weighting=True, model_directory=model_directory).load_ES_models()
    alpha_table = pd.read_csv(os.path.join(data_directory, 'alpha_table_m1.csv'), index_col=0)

    p = c.params()
    p.set_param("defaultKm", 0.00001) # M
    p.set_param("defaultVmax", 10) #mmol/gDw/hr
    p.set_param("maxCycles", 30)
    # p.set_param("maxCycles", 150)
    p.set_param("timeStep", 1)
    p.set_param('writeFluxLog', True)

    from scarcc.sim_engine.simulation_workflow import run_sim_workflow

    SG_list = ['folA', 'folP'] # Optional, Example: ['folA', 'folP', 'folC']
    DG_list = [['folA', 'folP']] # Optional, Example: [['folA', 'folP'], ['folC', 'folP'], ['folA', 'folC']]
    method_list = ['m1']
    simulation_kwargs = {
        'E0': E0,
        'S0': S0,
        'base': sim_chamber_directory,
        'p': p}


    df_container = run_sim_workflow(method_list=method_list, data_directory=data_directory,
                                    DG_list=DG_list, **simulation_kwargs)

if __name__ == '__main__':
    main()