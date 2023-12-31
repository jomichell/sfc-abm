import copy

class SimVars(dict):
    """Simple utility class to allow 'dot' notation for 
    accessing variables"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    
    def set(self, **kwargs):
        for var_name in kwargs:
            self[var_name] = kwargs[var_name]
        
    def __deepcopy__(self, memo=None):
        return SimVars(copy.deepcopy(dict(self), memo=memo))

    def subset(self, keys):
        return SimVars((k,self[k]) for k in keys if k in self)
    
    
class SimBloc:
    """A 'bloc' of a simulation model. This could be e.g. a 
    sector (firms, households) or an overall 'main model' bloc"""

    # Each 'bloc' has instance variables (i.e. per household)
    # and class variables (i.e. aggregated sector level variables)

    @classmethod
    def aggregate(cls, instance_list, var_list, lag=0):
        agg_vars = {}
        for var in var_list:
            agg_vars[var] = 0
            for instance in instance_list:
                agg_vars[var]+=instance.svars[lag][var]
        return agg_vars

    # allow direct 'dot' access on the class to the current
    # period state variables (decided against this,
    # better to use container-like access).

    # def __getattr__(self, item):
    #    return self.svars[0][item]

    # allow container-like access to the state variable dicts
    def __getitem__(self, idx):
        return self.svars[idx]
        
    def to_dict(self):
        return dict(self.svars[0])
    
    def __init__(self, model = None):

        # associate a bloc with a parent model if specified
        # and copy initial state variables. Parameters
        # are shared with parent models
        if model is not None:
            self.model = model
            self.svars = copy.deepcopy(model.svars)
            self.params = model.params
        else:
            # dict of state variables, indexed by lag
            self.svars = {}
            self.svars[0] = SimVars()
            self.svars[-1] = SimVars()
            self.params = SimVars()
            self.model = self

    def unpack(self):
        return (self.params,
                self.svars[0], self.svars[-1],
                self.model.bank)
        
    def set_svars(self, lag=0, **kwargs):
        self.svars[lag].set(**kwargs)

    def set_params(self, **kwargs):
        self.params.set(**kwargs)

    def get_svars(self, lag=0):
        return self.svars[lag]
        
    def lag_svars(self):
        self.svars[-1] = self.svars[0]
        self.svars[0] = SimVars()

    def incr(self, var, amount):
        if var in self.svars[0].keys():
            self.svars[0][var] += amount
        else:
            self.svars[0][var] = amount


