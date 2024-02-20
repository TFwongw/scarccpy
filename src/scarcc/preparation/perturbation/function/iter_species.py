"""iterate items to pass into function""" 

def iter_species(models,f,*args,**kwargs): 
    """Factory for passing extra_objects iterate with model into functions.
    
    Parameters
    ----------
    models : list or zip
        List of models or zip of models and extra_objects that will iterate together.
    f : function
        Function to be called with in put as model and extra_objects.
    *args : positional arguments following models or models and extra_objects
    **kwargs : keyword arguments

    Returns
    -------
    r_object : list
        List of return objects from function f.

    Examples
    --------
    # print model names
    >>> _ = iter_species([E0, S0], print)
    E0
    S0 
    # print model names and corresponding alpha as extra_objects
    >>> _ = iter_species(zip([E0, S0],[20,10]), print)
    E0 [20]
    S0 [10]
    """
    def simple_iter():
        for model in models:
            r_object.append(f(model,*args,**kwargs))

    # extra_object iterative with models
    def coupled_iter():
        for model, *extra_objects in models: 
            r_object.append(f(model,extra_objects,*args,**kwargs))

    r_object = []
    iter_fun = coupled_iter if isinstance(models, zip) else simple_iter
    iter_fun()
    return(r_object)
