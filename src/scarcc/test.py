# from .perturb.iter_species import iter_species
import sys
import os

# python_path = r'C:\Users\wongt\OneDrive\COMETS_Jupyter\Run_Comets\scarcc\src'
# sys.path.append(python_path)
# os.environ['PYTHONPATH'] = python_path

from scarcc import sim_preparation

# from scarcc.preparation.

from scarcc.perturb.iter_species import iter_species
##
from scarcc.fba_preparation import medium
from scarcc.sim_engine import result_processing

from scarcc.fba_preparation import basic_model, component
# from 
# from .test import test_basic_model

# print(pd.read_csv('../Data/_alpha_table_m3.csv').columns)


# from sim_engine import result_processing
