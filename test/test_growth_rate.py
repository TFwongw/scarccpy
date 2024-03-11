from scarcc.data_analysis import MethodDataFiller
import pandas as pd

df_sg = pd.read_csv('../Data/BM_SG_m1.csv', index_col=0)
df_dg = pd.read_csv('../Data/BM_DG_m1.csv', index_col=0)
df_flux = pd.read_csv('../Data/flux_analysis_m1.csv', index_col=0)
df_container = {'m1': {'SG': {'biomass': df_sg, 'flux': df_flux},
                        'DG': {'biomass': df_dg, 'flux': df_flux}}}
mdf = MethodDataFiller(df_container, data_directory='../Data')
mdf.fill_container()
mdf.write_to_csv()
