from dataclasses import dataclass

from typing import List, Union, Tuple, Dict, Any
import cometspy as c

from functools import partial
from scarcc.util import convert_arg_to_list
from scarcc.preparation import iter_species

@dataclass(kw_only=True)
class LayoutConfig:
    # medium related
    E0: "cobra.Model"
    S0: "cobra.Model"
    E_model: 'comets.model' = None
    S_model: 'comets.model' = None
    carbon_source_val: float = .1
    base_nutrients: List = None
    nutrients_val: int = 100
    
    # whether to generate layouts for corresponding culture
    co: bool = True
    mono: bool = True
    mono_E: bool = True
    mono_S: bool = True
    Smono_carbon_source: str = 'bulk_ac_e' 
    carbon_source_in_cultures: Dict[str, str] = None

    # comets model specifications
    initial_pop: float = 1e-8
    obj_style: str = 'MAX_OBJECTIVE_MIN_TOTAL'

    def __post_init__(self):
        if self.mono is False:
            self.mono_E, self.mono_S = False, False 
        self.Smono_carbon_source = 'bulk_ac_e' if any(['EX_bulk_ac_e' in ele.id for ele in self.S0.reactions]) else 'ac_e'
        if self.base_nutrients is None:
            self.base_nutrients = [
                "ca2_e", "cl_e", "cobalt2_e", "cu2_e","fe2_e", "fe3_e", "k_e", "mg2_e",
                "mn2_e", "mobd_e", "ni2_e", "o2_e", "pi_e", "so4_e", "zn2_e", "nh4_e"]
        if self.carbon_source_in_cultures is None:
            self.carbon_source_in_cultures = {'co': 'lcts_e', 'mono_E': 'lcts_e', 'mono_S': self.Smono_carbon_source}
        if E_model is None or S_model is None:
            self.set_comets_model()
    
    def set_comets_model(self):
        def create_comets_model(self, model):
            comets_model = c.model(model)
            comets_model.open_exchanges()
            comets_model.initial_pop = [0, 0, self.initial_pop]
            comets_model.obj_style = self.obj_style
            return comets_model
        self.E_model, self.S_model = iter_species([self.E0, self.S0], create_comets_model)

    def create_common_media(self, species: List['cobra.Model'], carbon_source: str, additonal_nutrients: List[str] = [''], additonal_nutrients_val: List[int] = [100]):
        additonal_nutrients, additonal_nutrients_val = convert_arg_to_list(additonal_nutrients), convert_arg_to_list(additonal_nutrients_val)
        l = c.layout(species)
        for nutrient in self.base_nutrients:
            l.set_specific_metabolite(nutrient, self.nutrients_val)
        if (additonal_nutrients != ['']):
            if (len(additonal_nutrients) == len(additonal_nutrients_val)):
                for add_nut, add_nut_val in zip(additonal_nutrients, additonal_nutrients_val): 
                    l.set_specific_metabolite(add_nut, add_nut_val)
            else:
                print(f'Set all additional nutrients to {additonal_nutrients_val[0]}')
                for add_nut in additonal_nutrients:
                    l.set_specific_metabolite(add_nut, additonal_nutrients_val[0])            

        l.set_specific_metabolite(carbon_source, self.carbon_source_val)
        return l

    def create_layout_object(self):
        partial_create_common_media = partial(self.create_common_media, carbon_source_val=self.carbon_source_val)
        co_layout = partial_create_common_media([self.E_model, self.S_model], carbon_source=self.carbon_source_in_cultures['co']) if self.co else None
        E0_layout = partial_create_common_media([self.E_model], carbon_source=self.carbon_source_in_cultures['mono_E'], 
                                                additonal_nutrients='met__L_e', additonal_nutrients_val=[100]) if self.mono_E else None
        # TODO: flux weight parameter ->10 
        Scarbon_source_val = self.carbon_source_val/10 if self.Smono_carbon_source=='bulk_ac_e' else self.carbon_source_val # ac_e as bulk 
        S0_layout = self.create_common_media([self.S_model], carbon_source=self.carbon_source_in_cultures['mono_S'], carbon_source_val=Scarbon_source_val) if self.mono_S else None
        return [co_layout, E0_layout, S0_layout]


@dataclass(kw_only=True)
class SimulationConfig(LayoutConfig):   
    current_gene: str
    p: 'comets.p'
    n_dir: str
    alpha_table: str
    checker_suffix: str = None
    return_sim: bool = False
    ko: bool = False
