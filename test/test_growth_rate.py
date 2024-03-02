from scarcc.data_analysis import MethodDataFiller
import pandas as pd

df_bm = pd.DataFrame({'cycle':[1,2,3], 'E0_folA_coculture': [4,5,6], 'E0_folP_coculture': [4,5,6]}).set_index('cycle')
df_flux = pd.DataFrame.from_dict({
    'folA': {'cycle': 1,
                'r1': 2,
                'r2': 3},
    'folP': {'cycle': 1,
                'r1': 2,
                'r2': 3}}, orient='index')


# mdf.fill_container()

df_sg = pd.read_csv('../Data/BM_SG_m1.csv', index_col=0)
df_dg = pd.read_csv('../Data/BM_DG_m1.csv', index_col=0)
df_flux = pd.read_csv('../Data/flux_analysis_m1.csv', index_col=0)
df_container = {'m1': {'SG': {'biomass': df_sg, 'flux': df_flux},
                        'DG': {'biomass': df_dg, 'flux': df_flux}}}
mdf = MethodDataFiller(df_container, data_directory='../Data')
mdf.fill_container()
mdf.write_to_csv()
