"""iterate items to passinto function"""
def iter_species(models,f,*args,**kwargs): 
    """Factory for passing extra_objects iterate with model into functions."""
    def simple_iter():
        for model in models:
            r_object.append(f(model,*args,**kwargs))

    # extra_object iterative with models
    def coupled_iter():
        for model, *extra_objects in models: 
            r_object.append(f(model,extra_objects,*args,**kwargs))

    r_object = list()
    iter_fun = coupled_iter if isinstance(models, zip) else simple_iter
    iter_fun()

    return(r_object)
