def create_common_media(Species, carbon_source = "lcts_e", carbon_source_val = 5e-3, 
                        nutrients_val = 100, add_nutrient = '', add_nutrient_val = [100]):    # carbon_source = 'lcts_e', 'ac_e' or 'glc__D_e'
    # add_nutrient = 'met__L_e' for E.coli monoculture
    
    # convert into list for enumeration 
    if type(add_nutrient)  == str:
        add_nutrient = [add_nutrient]
    if type(add_nutrient_val)  == int:
        add_nutrient_val = [add_nutrient_val]
        
    l = c.layout(Species)
    
    base_nutrients = ["ca2_e", "cl_e", "cobalt2_e", "cu2_e","fe2_e", "fe3_e", "k_e","mg2_e",
                  "mn2_e", "mobd_e", "ni2_e", "o2_e", "pi_e", "so4_e", "zn2_e","nh4_e"]
    for nutrient in base_nutrients:
        l.set_specific_metabolite(nutrient, nutrients_val)
        
    if (add_nutrient != ['']):
        if (len(add_nutrient) == len(add_nutrient_val)):
            for _,i in enumerate(zip(add_nutrient, add_nutrient_val)): 
                l.set_specific_metabolite(i[0], i[1])
        else:
            print(f'Set all additional nutrients to {add_nutrient_val[0]}')
            for _,i in enumerate(add_nutrient): 
                l.set_specific_metabolite(i, add_nutrient_val[0])
    
    l.set_specific_metabolite(carbon_source, carbon_source_val)  
    return(l)
    
def sim_cultures(model, layout, p=None, base = None, genes=''): # non_functional argument model
    # separate function into mono & coculture to prevent using wrong layer
    if type(layout) is list:
        layout = layout[0] # layout object stored inside list of one element unpacking from iter_species
    
    sim = c.comets(layout, p) 
    sim.working_dir = base
    print(sim.working_dir)
    
    try:
        sim.run()
        # if modify_bool==True:
        #     sim.total_biomass.rename(columns={'S0.ac':'S0.glc'},inplace=True)
    except:
        logging.exception(f"{sim.run_output}")
        print(f"{sim.run_output}")
    biomass_df = sim.total_biomass.set_index('cycle')
    biomass_df.columns =  rename_columns(biomass_df)
    return biomass_df, sim

def create_c(model, initial_pop=1e-8, obj_style='MAX_OBJECTIVE_MIN_TOTAL'):
    model = c.model(model)
    model.open_exchanges()
    model.initial_pop = [0, 0, initial_pop]
    model.obj_style = obj_style
    return(model)


def create_layout_object(E_model, S_model, carbon_source_val=5e-3, add_nutrient_val=[100], co_met=[0], co=True, mono=True, mono_S=True, Smono_carbon_source='bulk_ac_e'):
    add_nutrient_val = convert_arg_to_list(add_nutrient_val)
    partial_create_common_media = partial(create_common_media, carbon_source_val=carbon_source_val)
    
    # co_layout = partial_create_common_media([E_model, S_model], carbon_source='lcts_e') if co else None
    co_layout = partial_create_common_media([E_model, S_model], carbon_source='lcts_e', add_nutrient='met__L_e', add_nutrient_val=co_met) if co else None
    E0_layout = partial_create_common_media([E_model], carbon_source='lcts_e', add_nutrient='met__L_e', add_nutrient_val=add_nutrient_val) if mono else None
    
    Scarbon_source_val = carbon_source_val/10 if Smono_carbon_source=='bulk_ac_e' else carbon_source_val # ac_e as bulk 
    print('set Scarbon source as', Smono_carbon_source, Scarbon_source_val)
    
    S0_layout = create_common_media([S_model], carbon_source=Smono_carbon_source, carbon_source_val=Scarbon_source_val) if mono_S else None
    # return [co_layout, E0_layout, S0_ac_layout, S0_glc_layout]
    return [co_layout, E0_layout, S0_layout]