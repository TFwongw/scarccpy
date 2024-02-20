import os
def find_directory(desired_directory_name, script_filepath): # start with where the script is located
    """Find directory in the current directory or up maximum 2 levels of directory. 
    If not found, make directory at script location.
    
    Parameters
    ----------
    
    desired_directory_name : str
        Name of the directory to find or make.

    script_filepath : str
        Path to the script that is calling this function.

    Returns
    -------
    directory_path : str

    """

    def check_directory(directory_path):
        return os.path.isdir(directory_path)
    
    start_directory = os.path.abspath(script_filepath)
    current_directory = start_directory

    max_depth = 2
    while os.path.basename(current_directory) != 'scarccpy' and max_depth > 0:
        # Move up one level
        current_directory = os.path.dirname(current_directory)
        max_depth -= 1

        directory_path = os.path.join(current_directory, desired_directory_name)
        if check_directory(directory_path):
            return directory_path

    print(f'Directory {desired_directory_name} not exist, make directory at script level.')
    directory_path = os.path.join(os.path.dirname(start_directory), desired_directory_name)
    os.makedirs(directory_path)
    return directory_path
