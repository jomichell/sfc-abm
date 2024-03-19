from absteindl import Steindl
#from plotsims import plot
import logging

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# ----------------------------------------------------------
# Moderate stochastic/concentration balance, zeta = 0.77


# create the simulation model
sim1 = Steindl(num_firms = 1000, num_periods = 500, seed = 2)


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
    llambda = 0.75, # profit retention ratio (lambda reserved)
    r_L_bar = 0.03, # interest rate
    tau     = 0.25, # mark-up
    pr      = 1,    # labour productivity
    zeta    = 0.77   # size to revenue feedback 
)

# set the required initial variables
sim1.set_ivars(
    Y    = 34,
    Y_hr = 27,
    K    = 100,
    I    = 10,
    IV   = 5,
    F_n  = 2.8,
    F_r  = 2.1,
    r    = 0.06,
    u    = 1.46,
    D_h  = 80,
    D_f  = 136,
    L    = 80 + 136,
)

# initialise the model (creates firms with stochastic variation etc.)
sim1.initialise()

# run the model
sim1.run()

# generate plots
sim1.plot()


# to inspect results manually, uncomment:
# sim1.results

# -------------------------------------------------------------
# Strong stochastic effect, zeta = 0.5

# create a new sim based on sim1
sim2 = Steindl(num_firms = 1000, num_periods = 500, seed = 1)
sim2.copy_init(sim1)

# change a parameter
sim2.params.zeta = 0.5

sim2.initialise()
sim2.run()
sim2.plot()


# -------------------------------------------------------------
# Medium stochastic effect, zeta = 0.7

sim3 = Steindl(num_firms = 1000, num_periods = 500, seed = 2)
sim3.copy_init(sim1)

# change a parameter
sim3.params.zeta = 0.7

sim3.initialise()
sim3.run()
sim3.plot()


# -------------------------------------------------------------
# Limited stochastic effect, strong concentration, zeta = 0.9

sim4 = Steindl(num_firms = 1000, num_periods = 500, seed = 2)
sim4.copy_init(sim1)

# change a parameter
sim4.params.zeta = 0.9

sim4.initialise()
sim4.run()
sim4.plot()

 


# -----------------------------------------------------------
# Fully stochastic version, no concentration tendency


sim5 = Steindl(num_firms = 1000, num_periods = 500, seed = 1)
sim5.copy_init(sim1)

sim5.params.zeta = 0

sim5.initialise()
sim5.run()
sim5.plot()


# -----------------------------------------------------------
# 100% retained earnings, no distributed earnings

sim6 = Steindl(num_firms = 1000, num_periods = 500, seed = 1)
sim6.copy_init(sim1)

sim6.params.llambda = 1

sim6.initialise()
sim6.run()
sim6.plot()


# -----------------------------------------------------------
# Baseline model

sim7 = Steindl(num_firms = 1000, num_periods = 500, seed = 1)
sim7.copy_init(sim1)

# run a simplified 'baseline' model
sim7.params.baseline_model = True

sim7.params.llambda = 1
sim7.params.zeta = 1
#sim7.params.alpha2 = 0.2

sim7.initialise()
sim7.run()
sim7.plot()
