import os
import cometspy as c
import logging

from scarcc.preparation.metabolic_model import BasicModel
from scarcc.preparation.find_directory import find_directory
from scarcc.preparation.target_gene.gene_format_handler import get_DG_list, get_SG_list
from scarcc.preparation.perturbation.alpha_finder.monoculture import get_alpha_biomass_df
from scarcc.sim_engine.simulation_workflow import run_sim_workflow

logging.basicConfig(filename='scarcc_simulation.log',
                    level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')

def main():
    n_combos = None

    model_directory = find_directory('models', os.path.abspath(''))
    data_directory = find_directory('Data', os.path.abspath(''))
    sim_chamber_directory = find_directory('SimChamber', os.path.abspath(''))

    # initialize model
    E0, S0, all_components = BasicModel(flux_weighting=True,
                                    ac_scale=10, # Optional
                                    gal_scale=1/3, # Optional
                                    model_directory=model_directory).load_ES_models()
    DG_list = get_DG_list(os.path.join(data_directory, 'GeneCombos.csv'), n_combos=n_combos)
    SG_list = get_SG_list(DG_list)

    maf_kwargs = {
                'target_normalized_biomass': 0.5,
                'potential_genes': SG_list,
                'precision': 2,
                'acceptance_threshold_lower': 1}
    # default detailed_alpha_table=False
    alpha_table = get_alpha_biomass_df(model_list = [E0, S0], data_directory=data_directory, **maf_kwargs)

    p = c.params()
    p.set_param("defaultKm", 0.00001) # M
    p.set_param("defaultVmax", 10) #mmol/gDw/hr
    p.set_param("maxCycles", 550)
    # p.set_param("maxCycles", 30)
    p.set_param("timeStep", 1)
    p.set_param('writeFluxLog', True)

    method_list = ['m1']
    
    simulation_kwargs = {
        'E0': E0,
        'S0': S0,
        'base': sim_chamber_directory,
        'p': p,
        'mono': True,
        # 'save_raw_flux': True,
    }

    df_container = run_sim_workflow(method_list=method_list, data_directory=data_directory,
                                    SG_list=SG_list, DG_list=DG_list, **simulation_kwargs)

if __name__ == '__main__':
    main()