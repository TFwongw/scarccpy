# import
# import scarcc
from importlib import reload
import scarcc.preparation.find_directory
reload(scarcc.preparation.find_directory)
from scarcc.preparation.find_directory import find_directory

data_directory = find_directory('models', __file__)
print('out', data_directory)