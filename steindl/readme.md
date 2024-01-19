# ABSim #

A collection of Python tools for developing agent-based models

This is draft documention. This is a work in progress.


## Quickstart ##

Some quick and rough instructions on the main features of the model. For example simulation results, see: [simulations](./Steindl_ABM_simulations.md).

To create an agent such as firm with given parameters and state variables, we subclass `SimBloc`(see below):


```python
from absim import SimBloc

firm1 = SimBloc()

firm1.set_svars(y = 100, wb = 80)
firm1.set_svars(y = 75, wb = 40, lag = 1)
firm1.set_params(alpha = 0.1, gamma = 0.8)
firm1

```




    state vars:{0: {'y': 100, 'wb': 80}, -1: {'y': 75, 'wb': 40}}, params:{'alpha': 0.1, 'gamma': 0.8}



We can define specific features of agents by creating subclassing `SimBlock` and adding custom functions. These functions can take advantage of the `unpack` utility function to gain easy access to parameters 


```python
class Firm(SimBloc):
    def get_profit(self):
        c = self.get_svars()

        return c.y - c.wb


firm2 = Firm()
firm2.set_svars(y = 100, wb = 80)

firm2.get_profit()
```




    20



## In more detail ##

The module `absim` defines classes which provide basic functionality to simulate agent-based models.

### Model blocs ###


AB simulation models are thought of as being constructed from *blocs*, which are provided by the class `SimBloc`. Blocs are very loosely defined: a bloc can be used to represent a model, a sector or an agent. A bloc provides a way to collect and access *state variables* and *parameters*, and to define structure such as sectoral structure.

AB simulations are internally represented as a set of state variables which evolve through time. At any point, these only include those variables needed to continue simulation: past values of state variables are stored as 'results'. State variables can either be current or lagged variable values. For example, for a firm, previous period revenue might be a state variable because it affects current period investment.

Models always have a top-level bloc, which acts as a 'container' for the model as a whole, and provides functions to initialise the model, run simulations and other top-level functionality. Models may also contain lower-level blocs which can be used to collect together groups of agents of . For example, a two-sector model made up of households, firms may use a one 'bloc' to represent each household or firm plus an extra bloc for the overall model. Each agent will then have an independent set of state variables. These could be balance sheet items or 
current of lagged flow variables.



### Variable storage and access ###

The class `SimVars` is designed to overcome the problem of cumbersome algebraic syntax in models. In general, variables need to be kept out of the root namespace. But accessing variables in dictionaries is syntactically cumbersome and makes algebra difficult to read. For example, consider the macroeconomic identity $$Y = C + I$$ Ideally we'd like to be able to code this as `Y = C + I`, but this requires the three variables to defined in the root namespace. If the variables are stored in a dictionary, d, this becomes `d['Y'] = d['C'] + c['I']`. This is significantly more annoying to type, and also substantially harder to read. `SimVars` offers a compromise by allowing 'dot notation': 



```python
from absim import SimVars
d = SimVars()
d.C = 80
d.I = 20

d.Y = d.C + d.I

d.Y
```




    100



The class `SimBloc` sets up storage for state variables and model parameters using the `SimVar` class in the background. By default `SimBloc` assumes that current and one-period lagged variables constitute state variables.

We can create an object representing, for example, a firm, as follows


```python
from absim import SimBloc

firm = SimBloc()
firm.svars
```




    {0: {}, -1: {}}




```python

```
