
def get_BM_df(current_gene, n_dir, alpha_table, co=True, mono=True, mono_S=True, p=None, checker_suffix=None, E0=None, S0=None,return_sim=False, 
              ko=False, carbon_source_val=5e-3, add_nutrient_val=[100], initial_pop=1e-8, obj_style='MAX_OBJECTIVE_MIN_TOTAL', no_monoculture=2,
              reduce_SG_cycle=False, Smono_carbon_source='bulk_ac_e'):
                  
    # correction_factor = {'ac': ac_scale, 'gal': gal_scale, 'co2': co2_scale}
    # correction_factor = {k:v for k,v in correction_factor.items() if v is not None} if w else {}
    if p is None:
        sys.exit('Need to provide COMETS params')
    
    def get_p_short(p):
        p_short = copy.deepcopy(p)
        if p.get_param("maxCycles") > 200:
            p_short.set_param("maxCycles", 200) 
        return p_short
    
    genes=convert_arg_to_list(current_gene) 
    base = f"/panfs/jay/groups/0/harcombe/wong0755/comets_RPS/rep_{n_dir}/"
    # base = f"../comets_RPS/rep_{n_dir}/"
    if not os.path.exists(base):
        os.mkdir(base)  
        
    with E0 as m_E0, S0 as m_S0:
        M0 = [m_E0, m_S0] 
        if not ('Normal' in genes):
            alphas = iter_species(M0, get_alphas_from_tab, genes=genes, alpha_table=alpha_table)
            zip_arg = zip(M0, alphas)
            iter_species(zip_arg, alter_Sij,genes=genes, ko=ko)
        
        E_model,S_model = iter_species(M0, create_c, initial_pop=initial_pop, obj_style=obj_style) # create comet object with altered stiochiometry
        
        co_layout, *mono_layout = create_layout_object(E_model, S_model, carbon_source_val=carbon_source_val, add_nutrient_val=add_nutrient_val,
                                                        co=co, mono=mono, mono_S=mono_S, Smono_carbon_source=Smono_carbon_source)

        if co:
            print(p.all_params['maxCycles'], 'co_p')
            co_df, co_sim = sim_cultures([E_model, S_model], co_layout, base=base, p=p)
            full_df = co_df.add_suffix(f'_{genes}_coculture')
            
        else: 
            co_df, co_sim = None, None
            full_df = pd.DataFrame()
        
        if mono: 
            if mono_S is False:
                no_monoculture = 1
            zip_arg = zip([E_model, S_model], mono_layout[:no_monoculture]) # only E_model layout if mono_S is False
            
            # change maxCycle of param for monoculture to 200
            p_short = get_p_short(p) if ((len(current_gene) == 1) and reduce_SG_cycle) else p # only single gene inhibition 
            mono_output = iter_species(zip_arg, sim_cultures, base=base, p=p_short)
            mono_df, mono_sim = unpack_output_object(mono_output) 
            for new_df in mono_df:
                full_df = pd.concat([full_df, new_df.add_suffix(f'_{genes}_monoculture')],axis=1)
        else:
            mono_sim = None
        # out_dict = extract_dfs(mono_sim, co_sim, flux_correction_factor=correction_factor)
        out_dict = extract_dfs(mono_sim, co_sim, flux_correction_factor={})
        genes_str = '.'.join(genes)
        
        # adjust for checker board
        if checker_suffix:
            full_df = full_df.add_suffix(checker_suffix) # ((G1_lv, G2_lv)
            genes_str = genes_str + checker_suffix
        
        out_dict.update( {'Gene_inhibition': genes_str}) # for DG
        print(genes_str)
        
        full_df.columns = rename_columns(full_df) 

    return (full_df, out_dict) if not return_sim else (full_df, out_dict, co_sim)
