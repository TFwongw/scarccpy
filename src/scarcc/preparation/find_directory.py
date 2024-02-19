import os
def find_directory(desired_directory_name, script_filepath): # start with where the script is located
    """Find directory in the current directory or in the parent directory. 
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
        Path to the desired directory.
    """
    start_directory = os.path.dirname(script_filepath)
    print(os.path.dirname(start_directory))
    directory_path = os.path.join(start_directory, desired_directory_name)
    if not os.path.isdir(directory_path):
        print(f'Directory {desired_directory_name} not exist, make directory at script location')
        directory_path = os.path.join(os.path.abspath(start_directory), desired_directory_name) 
        os.makedirs(directory_path)
    return directory_path
