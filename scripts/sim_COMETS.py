import logging
import pandas as pd
import logging
import os
import pandas as pd
import itertools

import cometspy as c
from scarcc.preparation.find_directory import find_directory 
from scarcc.preparation.metabolic_model import BasicModel
from scarcc.sim_engine.simulation_workflow import CombinedAntibioticsSimulation

p = c.params()
p.set_param("defaultKm", 0.00001) # M 
p.set_param("defaultVmax", 10) #mmol/gDw/hr
p.set_param("maxCycles", 30)
# p.set_param("maxCycles", 150)
p.set_param("timeStep", 1)
p.set_param('writeFluxLog', True)

logging.basicConfig(filename='scarcc_simulation.log',
                    level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')

model_directory = find_directory('models', os.path.abspath(''))
data_directory = find_directory('Data', os.path.abspath(''))

chamber_directory = find_directory('SimChamber', os.path.abspath(''))

# initialize model
E0, S0, all_components = BasicModel(model_directory=model_directory, flux_weighting=True).load_ES_models()
alpha_table = pd.read_csv(os.path.join(data_directory, 'alpha_table_m1.csv'), index_col=0)

cas = CombinedAntibioticsSimulation(current_gene='folA', E0=E0, S0=S0, alpha_table=alpha_table, 
                                p=p, mono=True, base=chamber_directory)

biomass_df_list, flux_df_list = cas.get_BM_df()

pd.concat(biomass_df_list, axis=1).to_csv(os.path.join(data_directory, 'test_biomass_df.csv'))
pd.concat(itertools.chain(*flux_df_list)).to_csv(os.path.join(data_directory, 'test_flux_df.csv'))
