"""Create medium for the model for generating solution from FBA and MonocultureAlphaSearch"""

from collections.abc import Iterable
from typing import Union, List

def initialize_medium():
    """Returns a dictionary of default medium conditions for the model.
    Medium is without carbon source.
    10 as default value for COMETS maximum uptake rate"""
    nutrient_medium = {'EX_ca2_e': 10,
                        'EX_cl_e': 10,
                        'EX_cobalt2_e': 10,
                        'EX_cu2_e': 10,
                        'EX_fe2_e': 10,
                        'EX_fe3_e': 10,
                        'EX_k_e': 10,
                        'EX_mg2_e': 10,
                        'EX_mn2_e': 10,
                        'EX_mobd_e': 10,
                        'EX_ni2_e': 10,
                        'EX_o2_e': 10,
                        'EX_pi_e': 10,
                        'EX_so4_e': 10,
                        'EX_zn2_e': 10,
                        'EX_nh4_e': 10}
    return nutrient_medium

def change_medium(model: "cobra.Model", metabolite: Union["Metabolite" , List['Metabolite']], value: float, return_medium=False): # function for setting medium metabolite value 
    """Changes the value of a metabolite/ list of metabolite in the medium
    
    Parameters
    ----------
    model: cobra.Model
    metabolite: cobra.Metabolite
    value: float
    return_medium: bool
    
    Returns
    -------
    None or dict
    
    """
    medium = model.medium
    if not isinstance(value, Iterable):
        metabolite, value = [metabolite], [value]
    for m, v in zip(metabolite, value):
        medium[m] = v
    model.medium = medium
    if return_medium:
        return model.medium
