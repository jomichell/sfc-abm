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
    
    def __init__(self, model = None, lags = 1):
        self.lags = lags

        # associate a bloc with a parent model if specified
        # and copy initial state variables. Parameters
        # are shared with parent models
        if model is not None:
            self.model = model
            #self.svars = copy.deepcopy(model.svars)
            self.params = model.params
        else:
            self.params = SimVars()
            self.model = None
            
            # dict of state variables, indexed by lag
        self.svars = {}

        for lag in self.get_lags():
            self.svars[lag] = SimVars()

        # a bloc of variables which use for initialisation but
        # do not become state variables.
        self.ivars = SimVars()
        

    def __repr__(self):
        if self.model is None:
            return """state vars:{}, params:{}""".format(
                self.svars, self.params)
        else:
            return """state vars:{}""".format(
                self.svars)

    def get_lags(self):
        return list(range(0-self.lags, 1))
    
    def unpack(self):
        ''' return a tuple containing model parameters plus
        state vars for each lag '''

        lagged_svars = [self.svars[lag] for lag in self.get_lags()]
        lagged_svars.reverse()
        return tuple([self.params] + lagged_svars)
    
#                self.svars[0], self.svars[-1],
#                self.model.bank)
        
    def set_svars(self, lag=0, **kwargs):
        self.svars[0-lag].set(**kwargs)

    def set_params(self, **kwargs):
        self.params.set(**kwargs)

    def set_ivars(self, **kwargs):
        self.ivars.set(**kwargs)

        
    def get_svars(self, lag=0):
        return self.svars[lag]
        
    def lag_svars(self):
        #for lag in self.get_lags():
        for lag in self.get_lags()[:-1]:
            self.svars[lag] = self.svars[lag+1]

        self.svars[0] = SimVars()

    #def incr(self, var, amount):
    #    if var in self.svars[0].keys():
    #        self.svars[0][var] += amount
    #    else:
    #        self.svars[0][var] = amount


