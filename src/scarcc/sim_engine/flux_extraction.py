"""Extract biomass and flux data from cometspy solution object."""

import re
from dataclasses import dataclass
from typing import Dict, Union, List, Tuple
import itertools
import pandas as pd


import cobra
import cometspy as c

from scarcc.data_analysis.flux.flux_snapshot import adjust_flux_df
from scarcc.data_analysis import get_desired_cycle
from scarcc.utils import convert_arg_to_list

@dataclass
class SimObjectBase:
    """Base class for extracting biomass and flux data from cometspy solution object.
    
    Attributes
    ----------
    E0 : cobra.Model
    S0 : cobra.Model
    sim_object : c.solution
    alpha_table : pd.DataFrame
    current_gene : str    
    """

    E0: 'cobra.Model'
    S0: 'cobra.Model'
    sim_object: 'c.solution'
    alpha_table: pd.DataFrame
    current_gene: str

    def __post_init__(self):
        if isinstance(self.current_gene, list):
            self.current_gene = '.'.join(self.current_gene)
        self.biomass_df = (self.sim_object.total_biomass.set_index('cycle'))
        self.model_dict = {'E0': self.E0, 'S0': self.S0}
        self.species_list = [self.model_dict[ele] for ele in self.biomass_df.columns]
        self.culture = 'coculture' if len(self.species_list)>1 else f'monoculture'
        self.biomass_df = self.biomass_df.add_suffix(f'_{self.current_gene}_{self.culture}') # ?need to use rename_columns?

def get_flux_snapshot(sob: SimObjectBase, model: 'cobra.Model' = None): 
    """Get flux snapshot for a model in a simulation object."""
    model_id = model.id
    biomass_df = sob.biomass_df.filter(regex=f'{model_id}') # suffix already added
    flux_df = sob.sim_object.fluxes_by_species[f'{model_id}'].copy()
    flux_df['Species'] = model_id
    flux_df['culture'] = sob.culture
    flux_df['Gene_inhibition'] = sob.current_gene
    flux_df = adjust_flux_df(model, flux_df, sob.current_gene, alpha_table=sob.alpha_table)
    
    # get the last cycle
    desired_cycle = get_desired_cycle(biomass_df)
    
    # TODO: snapshot handler
    snap_shot_cycle = int(desired_cycle.cycle_max_gr.iloc[0])
    snap_shot = flux_df.query('cycle == @snap_shot_cycle')
    return snap_shot

def extract_biomass_flux_df(E0: 'cobra.Model', S0: 'cobra.Model', sim_object_list: List['c.solution'],
                            alpha_table: pd.DataFrame, current_gene: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Extract biomass and flux dataframes from a list of cometspy solution objects.
    
    Parameters
    ----------
    E0 : cobra.Model
    S0 : cobra.Model
    sim_object_list : List[c.solution]
        [solution_co, solution_mono_E, solution_mono_S]
    alpha_table : pd.DataFrame
        subset of alpha table(specific lv_pair in checkerboard) or full alpha table 
    current_gene : str

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        biomass_df, flux_df
        serve as output for the get_BM_df function
    """

    def extract_biomass_flux_df_per_sim(sob: SimObjectBase) -> Dict[str, Union[pd.DataFrame, List]]: 
        """Extract biomass and flux dataframes per simulation object.
        
        Parameters
        ----------
        sob : SimObjectBase
        
        Returns
        -------
        Dict[str, Union[pd.DataFrame, List]]
            will concat altogether altogether for efficiency
        """
        biomass_flux_dict = dict()
        biomass_flux_dict['biomass'] = sob.biomass_df
        # biomass_flux_dict['flux'] = pd.concat([get_flux_snapshot(sob, model_id=ele) for ele in sob.species_list])
        biomass_flux_dict['flux'] = [get_flux_snapshot(sob, model=ele) for ele in sob.species_list]
        return biomass_flux_dict['biomass'], biomass_flux_dict['flux']
    
    def construct_sob(sim_object: 'comets.solution'):
        """Construct SimObjectBase object."""
        return SimObjectBase(E0=E0, S0=S0, sim_object=sim_object, alpha_table=alpha_table, current_gene=current_gene)
    
    # extraction workflow
    biomass_list, flux_nested_list = zip(*[extract_biomass_flux_df_per_sim(construct_sob(sim_object)) for sim_object in convert_arg_to_list(sim_object_list)])
    biomass_df = pd.concat(biomass_list, axis=1)
    flux_df = pd.concat(itertools.chain(*flux_nested_list), axis=0) # concatenated for all co and mono sim_objects
    return biomass_df, flux_df

def rename_columns(df):
    """Rename columns of a data frame."""
    df.columns = [re.sub('S0_ac_','S0.ac_', ele) for ele in df] # S0_ac -> S0.ac
    df.columns = [re.sub('S0_gal_','S0.gal_', ele) for ele in df] # S0_ac -> S0.ac
    df.columns = [re.sub(',','.',
                    re.sub('\'|\(|\)| |\[|\]','',ele)) # ('gene1', 'gene2') -> gene1.gene2
                    for ele in df.columns]
    return(df.columns)
