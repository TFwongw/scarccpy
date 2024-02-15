from importlib import reload
import logging

logging.basicConfig(filename='test_basic_model.log', 
                    level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')

from scarcc.preparation.metabolic_model import BasicModel

E0, S0, all_components = BasicModel(flux_weighting=True).load_ES_models()
logging.debug('Acetate pp Reaction. %s', str(E0.reactions.ACt2rpp))
logging.debug('Acetate exchange Reaction. %s', str(E0.reactions.EX_bulk_ac_e))
logging.debug('Galactose producing Reaction. %s', str(E0.reactions.LACZ))
