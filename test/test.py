# from .perturb.iter_species import iter_species
# import sys
# import os

# # python_path = r'C:\Users\wongt\OneDrive\COMETS_Jupyter\Run_Comets\scarcc\src'
# # sys.path.append(python_path)
# # os.environ['PYTHONPATH'] = python_path

# from scarcc import sim_preparation

# # from scarcc.preparation.

# from scarcc.perturb.iter_species import iter_species
# ##
# from scarcc.fba_preparation import medium
# from scarcc.sim_engine import flux_extraction

# from scarcc.fba_preparation import basic_model, component
# # from 
# # from .test import test_basic_model

# # print(pd.read_csv('../Data/_alpha_table_m3.csv').columns)


# # from sim_engine import result_processing
# from scarcc.sim_engine.output import A
import pandas as pd

df_bm = pd.DataFrame({'cycle':[1,2,3], 'E0_folA_coculture': [4,5,6], 'E0_folP_coculture': [4,5,6]}).set_index('cycle')
df_flux = pd.DataFrame.from_dict({
    'folA': {'cycle': 1,
                'r1': 2,
                'r2': 3},
    'folP': {'cycle': 1,
                'r1': 2,
                'r2': 3}}, orient='index')

class A:
    def get_result():
        return df_bm, df_flux


def run_A(x):
    a = A()
    return A.get_result()

def concat_result(result_dict):
    def concat_df(key, df_list):
        if 'biomass' in key:
            return pd.concat(df_list, axis=1)
        return pd.concat(df_list, axis=0)
    return {k: concat_df(k, v) for k, v in result_dict.items()}

def unpack_future_result_per_key(result_list):
    keys = ['biomass', 'flux']
    # return {k: v for k, v in zip(keys, zip(*[r.result() for r in result_list]))}
    return dict(zip(keys, zip(*result_list)))

import concurrent.futures
if __name__ == '__main__':
    task_dict = dict()
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        for key in ['x']:
            task_dict[key] = [executor.submit(run_A, key) for i in range(2)]
    df_container = {key: [future.result() for future in future_list]
                    for key, future_list in task_dict.items()}
    df_container = {key: unpack_future_result_per_key(result_list) for key, result_list in df_container.items()}
    
    x = df_container['x']['biomass']
    y = df_container['x']['flux']
    # xx = pd.concat(x)
    # yy = pd.concat(y)
    # print(xx)
    # print(yy)
    df_container = {k: concat_result(sub_container) for k, sub_container in df_container.items()} # column-wise for biomass, row-wise for flux
    print(df_container['x']['biomass'])
    print(df_container['x']['flux'])

    print(df_container)