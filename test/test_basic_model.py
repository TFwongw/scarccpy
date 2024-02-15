# from scarcc.preparation.core import flux_weighting

# from scarcc.fba_preparation.basic_model import basic_model
# from scarcc.fba_preparation import component
# from scarcc.fba_preparation import medium

import logging
# configure logging
logging.basicConfig(filename='test.log', 
                    level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')


# from scarcc.preparation.metabolic_model import basic_model 
# from scarcc.preparation.metabolic_model.basic_model import basic_model 

import scarcc
from importlib import reload
reload(scarcc)
from scarcc.preparation.metabolic_model.basic_model import basic_model 
# from scarcc.preparation.metabolic_model.core import component
s = basic_model().load_ES_models()


