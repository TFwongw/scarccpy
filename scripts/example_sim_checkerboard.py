import os
import pandas as pd
import numpy as np  
import cometspy as c
import ast

from scarcc.preparation.metabolic_model import BasicModel
from scarcc.preparation.find_directory import find_directory 
from scarcc.sim_engine.checkerboard_workflow import run_checkerboard_workflow

if __name__ == "__main__":
    # get file directory 
    model_directory = find_directory('models', os.path.abspath(''))
    data_directory = find_directory('Data', os.path.abspath(''))
    sim_chamber_directory = find_directory('SimChamber', os.path.abspath(''))

    # initialize model
    E0, S0, all_components = BasicModel(model_directory=model_directory, flux_weighting=True).load_ES_models()

    # read alpha_table_checkerboard
    alpha_table_path = os.path.join(data_directory, 'alpha_table_checkerboard.csv')
    alpha_table_checkerboard = pd.read_csv(alpha_table_path, index_col=0, converters={'lv_pairs': ast.literal_eval})
    query_pairs = [(0, 0), (1,0), (0,2), (1,2)] # Uncomment this for full checkerboard
    alpha_table_checkerboard = alpha_table_checkerboard.query('lv_pairs in @query_pairs') # Uncomment this for full checkerboard

    p = c.params()
    p.set_param("defaultKm", 0.00001) # M
    p.set_param("defaultVmax", 10) #mmol/gDw/hr
    p.set_param("maxCycles", 30)
    # p.set_param("maxCycles", 550) # uncomment this for full simulation
    p.set_param("timeStep", 1)
    p.set_param('writeFluxLog', True)

    checkerboard_simulation_kwargs = {
        'E0': E0,
        'S0': S0,
        'base': sim_chamber_directory,
        'p': p}

    df_container = run_checkerboard_workflow(
        alpha_table = alpha_table_checkerboard,
        data_directory = data_directory,
        **checkerboard_simulation_kwargs)