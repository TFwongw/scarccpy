import re
from data_analysis.flux import flux_correction

def unpack_output_object(output):
    df, *sim_object = zip(*output)
    return df, sim_object
# df, sim_object = unpack_output_object(mono_output)

def extract_dfs_from_sim_object(sim_object, flux_correction_factor):
    species_name = [ele for ele in sim_object.total_biomass.columns[1:]]
    if len(species_name)>1:
        culture = 'coculture' 
        out_dict = {f'{culture}_media' : sim_object.media}
    else: 
        culture = f'monoculture'
        out_dict = {f'{species_name[0]}_{culture}_media' : sim_object.media}
    for species in species_name:
        flux_df =  sim_object.fluxes_by_species[f'{species}']
        if species == 'E0':
            print('correct E')
            flux_df = flux_correction(flux_correction_factor, flux_df)
            print('maxac', flux_df['EX_bulk_ac_e'].abs().max())
        out_dict[f'{species}_{culture}_flux'] = flux_df
    return out_dict
        
# unpack_output_object
def extract_dfs(mono_sims, co_sim, flux_correction_factor=None):
    out_dict = extract_dfs_from_sim_object(co_sim, flux_correction_factor) if co_sim else dict() # initialize out_dict from coculture if have co_sim object 
    if mono_sims:
        for sim_object in mono_sims[0]:
            out_dict.update(extract_dfs_from_sim_object(sim_object, flux_correction_factor))
    out_dict = {k: v.to_dict() for k,v in out_dict.items()}
    return out_dict

def rename_columns(df):
    df.columns = [re.sub('S0_ac_','S0.ac_', ele) for ele in df] # S0_ac -> S0.ac
    df.columns = [re.sub('S0_gal_','S0.gal_', ele) for ele in df] # S0_ac -> S0.ac
    df.columns = [re.sub(',','.',
           re.sub('\'|\(|\)| |\[|\]','',ele)) # ('gene1', 'gene2') -> gene1.gene2
           for ele in df.columns]
    return(df.columns)
    
def gene_index_culture_col_df(analysis_df): 
    analysis_df['Gene_inhibition'] =  ['.'.join(map(str, convert_arg_to_list(l))) for l in analysis_df.Gene_inhibition] # SG, DG, checkerboard g1.g2 form
    analysis_df = analysis_df.set_index('Gene_inhibition')
    return analysis_df
