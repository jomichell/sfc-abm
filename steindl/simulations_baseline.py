# Attempt to create a set of 'baseline' parameterisations and scenarios

from absteindl import Steindl
#from plotsims import plot
import logging

logger = logging.getLogger()
logger.setLevel(logging.INFO)


# create the simulation model
sim1 = Steindl(num_firms = 10, num_periods = 500)

sim1.flags.baseline_model = True
sim1.flags.firm_results = True

# set the model parameters
sim1.set_params(
    alpha1  = 0.7,  # consumption out of income
    alpha2  = 0.1,  # consumption out of wealth
    gamma0  = 0.02, # trend growth
    gamma_m = 0.05,  # profit margin sensitivity
    gamma_u = 0.05, # utilisation rate sensitivity
    iota    = 0.2,  # desired inventory to output ratio
    theta   = 0.5,  # desired liquidity ratio
    v       = 4,    # capital full output ratio
    llambda = 1, # profit retention ratio (lambda reserved)
    r_L_bar = 0.03, # interest rate
    tau_bar = 0.25, # mark-up
    kappa   = 3,   # degree of mark-up adjustment to market share
    pr      = 1,    # labour productivity
    zeta    = 0.98   # size to revenue feedback 
)

# set the required initial variables
sim1.set_ivars(
    Y    = 34,
    Y_h  = 27,
    K    = 100,
    I    = 10,
    IV   = 5,
    F_n  = 2.8,
    F_r  = 2.1,
    r    = 0.06,
    u    = 1.46,
    D_h  = 100,
    D_f  = 0,
    L    = 80,
)

# initialise the model (creates firms with stochastic variation etc.)
sim1.initialise()

# run the model
sim1.run()

# generate plots
sim1.plot(skip_periods = 50)


# create a new sim based on sim1
sim3 = Steindl(num_firms = 100, num_periods = 1000, seed = 1)
sim3.copy_init(sim1)

sim3.params.zeta = 0.98
sim3.baseline_model = False


sim3.initialise()
sim3.run()
sim3.plot()
