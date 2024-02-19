from scarcc.preparation.find_directory import find_directory

# data_directory = find_directory(os.path.abspath(__file__), 'models')
data_directory = find_directory('models', __file__)
data_directory