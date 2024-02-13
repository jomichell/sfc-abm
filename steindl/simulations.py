from absteindl import Steindl
#from plotsims import plot
import logging

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# ----------------------------------------------------------
# zeta = 0.77

sim1 = Steindl(num_firms = 1000, num_periods = 500, seed = 2)

sim1.set_params(
    alpha1  = 0.7,  # consumption out of income
    alpha2  = 0.1,  # consumption out of wealth
    alpha3  = 0.0,  # inelastic saving out of salaries etc.
    gamma0  = 0.02, # trend growth
    gamma_r = 0.2,  # profit rate multiplier
    gamma_u = 0.04, # utilisation rate multplier
    iota    = 0.2,  # desired inventory to output ratio
    theta   = 0.5,  # desired liquidity ratio
    v       = 4,    # capital full output ratio
    llambda = 0.75, # profit retention ratio (lambda reserved)
    r_l_bar = 0.03, # interest rate
    tau     = 0.25, # mark-up
    pr      = 1,    # labour productivity
    zeta    = 0.77   # size to revenue feedback 
)

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

sim1.initialise()
sim1.run()

sim1.plot()

# -------------------------------------------------------------
# zeta = 0.5

# create a new sim based on sim1
sim2 = Steindl(num_firms = 1000, num_periods = 500, seed = 1)
sim2.copy_init(sim1)

# change a parameter
sim2.params.zeta = 0.5

sim2.initialise()
sim2.run()
sim2.plot()


# -------------------------------------------------------------
# zeta = 0.7

sim3 = Steindl(num_firms = 1000, num_periods = 500, seed = 2)
sim3.copy_init(sim1)

# change a parameter
sim3.params.zeta = 0.7

sim3.initialise()
sim3.run()
sim3.plot()


# -------------------------------------------------------------
# zeta = 0.9

sim4 = Steindl(num_firms = 1000, num_periods = 500, seed = 2)
sim4.copy_init(sim1)

# change a parameter
sim4.params.zeta = 0.9

sim4.initialise()
sim4.run()
sim4.plot()

 
