# SFC ABM steindl model from this paper:
# https://www.postkeynesian.net/working-papers/1412/
#
# The code used to generate the simulations shown in the
# paper contained some accounting errors. These errors are
# corrected in this version of the code.

### Set working directory
import os
os.chdir("C:/Users\Karsten Kohler/Dropbox/ABM project with Jo")

# Load libraries
import matplotlib.pyplot as pyplot
import pandas as pd
import numpy as np
import sys, logging
import time

# Start the timer
start_time = time.time()

# Change the below to DEBUG for more logging output (creates a log file)
logging.basicConfig(level=logging.INFO)

# Set a seed
np.random.seed(42)

##############################################################################
### %% START INITIALISATION
#############################################################################

start_year = 1
end_year   = 8000
num_years  = end_year - start_year
years = range(start_year, end_year)

# Number of firms
num_firms = 50

# Parameters

# Household consumption propensities
alpha1 = 0.7 # consumption out of income
alpha2 = 0.1 # consumption out of wealth
alpha3 = 0.0 # inelastic saving out of salaries etc.

# Investment behaviour
gamma0 =  0.02 # trend growth
gamma_r = 0.2  # profit rate multiplier
gamma_u = 0.04  # utilisation rate multplier

iota    = 0.2  # desired inventory to output ratio
theta   = 0.5  # desired liquidity ratio
v       = 1 #4 # capital full output ratio

llambda = 1#0.75#1 # profit retention ratio
# NB lambda is a reserved word in python.

# Exogenous
r_l_bar = 0.03
tau     = 0.25 # mark-up
pr      = 1    # labour productivity

# strength of effect of capital
# share on aggregate demand allocation
zeta    = 1  # 0.7

logging.info("initialising")

# Set up dict to collect aggregate time series
ts = dict()
for var_name in ['c', 'k', 'l', 'f_t', 'f_r', 'f_d', 'f_n',
                 'g_i', 'g_y', 'i', 'iv', 'y_s', 'y', 'y_hr', 'y_f',
                 'm', 'm_h', 'm_f', 'm_t', 'r', 'u', 'wb', 'it',
                 'r_l', 'r_m']:
    ts[var_name] = pd.Series(index=years)

# extra series to store summary variables which are not used
# in the model
xts = dict()
for var_name in ['hedge', 'spec', 'ponzi', 'bankrupt']:
    xts[var_name] = pd.Series(index=years)

# create list of time-series 'balance sheets' for individual firms
firm_ts = []

for firm in range(num_firms):
    bs_series = {}               # create new dictionary for each firm
    for key in ts.keys():         # loop over entries in ts dictionary (called "keys")
        # Duplicate the aggregate balance sheet variables into the firm_ts list
        bs_series[key]      =  pd.Series(index=years)

    # plus some additional summary variables
    bs_series['k_share'] = pd.Series(index=years)
    bs_series['y_share'] = pd.Series(index=years)
    bs_series['size']  =   pd.Series(index=years)     
    firm_ts.append(bs_series)   # add balance sheet dictionary to firm_ts list 

### Set up initial values

# Required for macro solution
ts['y'][start_year] = 34.08
ts['y_hr'][start_year]  = 27.2
ts['k'][start_year]   = 100
ts['iv'][start_year]   = 5

ts['m_h'][start_year] = 79.5
ts['m_f'][start_year] = 132.6
ts['l'][start_year] = \
    ts['m_h'][start_year] + ts['m_f'][start_year]
ts['f_d'][start_year]   = 0.704

# Used for firm level initialisation
ts['i'][start_year]   = 10.28
ts['it'][start_year]   = 3.59
ts['wb'][start_year]    = 23.3 
ts['f_r'][start_year]   = 2.126
ts['f_n'][start_year]   = 2.83

ts['r'][start_year]    = 0.0583
ts['u'][start_year]    = 1.457

logging.debug("""initial variable values: 
y: {}
y_hr: {}
iv: {}
l: {}
m_f: {}
m_h: {}
""".format(
    ts['y'][start_year],
    ts['y_hr'][start_year],
    ts['iv'][start_year],
    ts['l'][start_year],
    ts['m_f'][start_year],
    ts['m_h'][start_year]))

logging.debug("allocating initial firms' balance sheets")

# Now set up the (sparse) firm balance sheets for the initial year with
# an initial capital allocation and duplicates of the aggregates for
# rate of profit and capacity utilisation (this will result in slight faulty
# investment demand decisions by firms in the first year but won't result in
# stock-flow inconsistency).

sy = start_year
total_k_rnd =0

for firm in firm_ts:  # iterate over each element in the firm_ts list
    firm['k_rnd'] = np.random.random()   # assign random numbers to firms to allocate capital
    total_k_rnd += firm['k_rnd']          # accumulate capital allocation into a running total

    # now we know the total of the random numbers, we can work out the relative size
    # of a firm and used that to distribute the aggregate initial values over the firms
for firm in firm_ts:  

    # work out the proportion of total capital this firm gets
    firm['k_share'][sy] = firm['k_rnd'] / total_k_rnd

    # and how this compares to the 'average' allocation, to give
    # the relative 'size' of the firm 
    firm['size'][sy] = firm['k_share'][sy] * num_firms

    # allocate an initial capital stock to the firm based on total initial capital stock k
    firm['k'][sy] =  ts['k'][sy] * firm['k_share'][sy]

    # and an initial share of demand
    firm['y'][sy] =  ts['y'][sy] * firm['k_share'][sy]

    # set u and r equal to the aggregate figures
    firm['u'][sy] =  ts['u'][sy]
    firm['r'][sy] =  ts['r'][sy]

    # allocate inventories, debt and cash on the basis of firm size
    firm['m_f'][sy]  =  ts['m_f'][sy] * firm['k_share'][sy]
    firm['l'][sy]    =  ts['l'][sy] * firm['k_share'][sy] 

    firm['it'][sy]   =  ts['it'][sy] * firm['k_share'][sy]
    firm['wb'][sy] =  ts['wb'][sy] * firm['k_share'][sy]

    firm['f_n'][sy] =  ts['f_n'][sy] * firm['k_share'][sy]
    firm['f_r'][sy] =  ts['f_r'][sy] * firm['k_share'][sy]    

    firm['i'][sy] =  ts['i'][sy] * firm['k_share'][sy]
    firm['iv'][sy] =   ts['iv'][sy] * firm['k_share'][sy]             

#############################################################
### %%  Main loop of iterative model solution
############################################################

logging.info("starting simulation for {}".format(years))

    # Main solution loop has three steps
    # 1) firm-level decision-making based on previous-year values.
    # 2) aggregate firm level decisions and use these to solve macro model.
    # 3) allocate macro magnitudes among firms according to market shares.

for year in years[1:]:

    logging.info("start of year {:d}".format(year))

    # some aggregate values are (weighted) totals of individual values
    agg_m        = 0
    agg_i        = 0
    agg_l        = 0
    agg_y_s      = 0 
    agg_wb       = 0 # wage bill
    agg_Dl_d     = 0 # change in loans
    agg_f_lb     = 0 # predicted current costs
    agg_lq_d     = 0 # desired minimum firms' deposits 
    agg_lq_e     = 0 # expected firm's deposits
    agg_iv_d     = 0 # desired inventories
    agg_f_t_e    = 0 # predicted profits
    
####%%%  1) Firm level decision making 
    logging.debug("*** Firm-level decision making ****")

    # A stochastic element is used to distribute demand across firms:
    # we assign a random 'share' magnitude to each firm which is then used
    # to distributing aggregate demand and other macro quantities

    total_rnd = 0

    for firm in firm_ts:
        # first assign random numbers to firms to create random element in the allocation of output 
        firm['rnd'] = np.random.random()
        total_rnd += firm['rnd']
    
    for firm in firm_ts:    
        # calculate a pure stochastic 'share' for each firm
        s_share = firm['rnd'] / total_rnd
        
        # use previous-year 'firm size' and 'capital share' to determine
        # mark-up and distribution of goods demand for current year
        fsize  = firm['size'][year-1]
        kshare = firm['k_share'][year-1]

        # use a linear combination of capital size and stochastic 'share'
        # to give the share of total expenditure that this firm will receive
        firm['y_share'][year] = (zeta * kshare) + ((1-zeta) * s_share) ## eq. 29

        # Investment (expenditure) function in growth terms
        firm['g_i'][year] = gamma0 + (gamma_r * firm['r'][year-1]) \
            + (gamma_u * firm['u'][year-1])  ## eq. 2

        # Investment (expenditure) in level terms
        firm['i'][year] = firm['g_i'][year] * firm['k'][year-1]

        # Sum over all investment expenditures to get aggregate figure
        agg_i += firm['i'][year]

        # We assume that larger firms can impose a greater mark-up
        f_tau  = tau * fsize ## eq 11

        # Thus a greater gross profit margin
        firm['m'][year] = f_tau/(1+f_tau)  ## eq. 12
        
        # Weighted average of margins gives us aggregate (is this right? KK: yes)
        agg_m += firm['y_share'][year] * firm['m'][year]        

        # Firm decides how much to produce on the basis of
        # previous year's sales, intended growth of
        # capital stock and desired inventories

        # Expected demand is previous year's demand
        # scaled by growth of capital stock
        
        y_e = (firm['y'][year-1] * (1+firm['g_i'][year])) ## eq. 7

        # desired inventories calculated on basis
        # of expected sales
        iv_d = iota * y_e ## eq. 8
        agg_iv_d += iv_d 

        # Production based on predicted sales and 
        # current excess inventories
        firm['y_s'][year] = iv_d - firm['iv'][year-1] + y_e ## eq. 9
        # rule out negative production
        if firm['y_s'][year] < 0:
            firm['y_s'][year] = 0
        
        agg_y_s += firm['y_s'][year]

        # Wage bill = value of production y_s less proft margin
        # times that value of production ## eq. 13
        firm['wb'][year]   = firm['y_s'][year] - \
            (firm['m'][year] * firm['y_s'][year])
        
        agg_wb+= firm['wb'][year]
        
        # Net interest payments
        # NB this assumes that firms already know the rate of interest in the
        # _current_ period. Will not be the case once banks behaviour is
        # modelled in more detail

        r_l_tmp = r_m_tmp = r_l_bar
        
        firm['it'][year] =  (
            r_l_tmp * firm['l'][year-1]) - (
            r_m_tmp * firm['m_f'][year-1])

        # Demand for loans is determined on the basis of desired
        # liquidity and intended investment

        # Predicted current costs
        # See note above on 'net interest payments'

        f_lb = (firm['wb'][year] + firm['it'][year])
      
        agg_f_lb += f_lb
        
        # Predicted net profits is predicted sales minus predicted costs
        f_n_e = y_e - f_lb

        # Predicted retained profits (not sure this is reported in the paper?)
        if f_n_e > 0:
            f_r_e = llambda * f_n_e
        else:
            f_r_e = f_n_e

        # Assume desired minimum liquidity is a proportion of current costs
        lq_d = theta * f_lb #  eq. 14
      
        agg_lq_d += lq_d

        # Calculate expected liquidity this period on the basis of
        # expected profits

        lq_e = firm['m_f'][year-1] + f_r_e - firm['i'][year]    ## eq. 15  
        
        agg_lq_e += lq_e

        # Compare expected to desired liquidity and
        # decide on loan demand accordingly

        if lq_e > lq_d: ## eq. 16
            firm['l'][year] = firm['l'][year-1]
        else:
            firm['l'][year] = firm['l'][year-1] + lq_d - lq_e

        # Sum loans outstanding to get total
        agg_l  += firm['l'][year]
        agg_Dl_d += firm['l'][year] - firm['l'][year-1]

        logging.log(8, """Firm -
        m_f(-1) : {:f}
        g_i     : {:f} 
        i     : {:f}
        y_e     : {:f}
        it(-1)  : {:f}
        f_lb    : {:f}
        f_n_e   : {:f}
        f_r_e   : {:f}
        lq_d    : {:f}
        lq_e    : {:f}
        Dl    : {:f}
        l     : {:f}      
        """.format(
            firm['m_f'][year-1],
            firm['g_i'][year],
            firm['i'][year],
            y_e,
            firm['it'][year-1],
            f_lb,
            f_n_e,
            f_r_e,
            lq_d,
            lq_e,
            firm['l'][year] - firm['l'][year-1],
            firm['l'][year]         
        ))

    logging.debug("""Inserting aggregates from micro:
    y_s: {}
    m  : {}
    i: {}
    wb: {} 
    l {}   
    """.format(
        agg_y_s,
        agg_m,
        agg_i,
        agg_wb,
        agg_l
        ))
    
    ## Save aggregated values in ts list
    ts['i'][year]   = agg_i
    ts['y_s'][year] = agg_y_s
    ts['wb'][year]  = agg_wb
    ts['m'][year]   = agg_m
    ts['l'][year]   = agg_l

    logging.debug(""" summary firm level aggregates:
    Dl_d : {:f}
    f_lb : {:f}
    lq_d : {:f}     
    lq_e : {:f}
    iv_d : {:f}
    """.format(
        agg_Dl_d,
        agg_f_lb,
        agg_lq_d,
        agg_lq_e,
        agg_iv_d
        ))
  
######### %%% 2) Solve the macro model to get aggregate numbers ###################
    logging.debug("solving aggregate model for year {}".format(year))
    
    r_l = r_l_bar
    r_m = r_l
    ts['c'][year] = alpha1 * ts['y_hr'][year-1] + alpha2 * ts['m_h'][year-1]  # cons spending
    ts['i'][year] = ts['i'][year]  # investment spending
    ts['y'][year] = ts['c'][year] + ts['i'][year]  # total expenditure
    ts['y_s'][year] = ts['y_s'][year]  # total production (supply side)
    ts['k'][year] = ts['k'][year-1] + ts['i'][year]  # production of capital goods (demand determined)
    ts['iv'][year] = ts['iv'][year-1] + (ts['y_s'][year] - ts['y'][year])  # chg in firms' inventories
    ts['it'][year] = (r_l * ts['l'][year-1]) - (r_m * ts['m_f'][year-1])  # net interest payments
    ts['wb'][year] = ts['wb'][year]  # Wage bill
    ts['f_t'][year] = ts['y'][year] - ts['wb'][year]  # gross profits
    ts['f_n'][year] = ts['f_t'][year] - ts['it'][year]  # net profits
    ts['l'][year] = ts['l'][year]  # firms sector loan demand
    
    # summary variables (not used for model solution)
    ts['g_y'][year] = (ts['y'][year] - ts['y'][year-1]) / ts['y'][year-1]  # GDP growth rate
    ts['g_i'][year] = (ts['i'][year] / ts['i'][year-1]) - 1  # investment growth rate 
    ts['r'][year] = ts['f_t'][year] / ts['k'][year]  # aggregate profit rate
    ts['y_f'][year] = ts['k'][year] / v  # aggregate capacity utilisation
    ts['u'][year] = ts['y'][year] / ts['y_f'][year]  # capacity utilisation
    ts['m'][year] = (ts['y'][year] - ts['wb'][year]) / ts['y'][year]  # aggregate profit share
    
    logging.debug("""macro solutions
    y    :{:f}
    y_s  :{:f}
    iv   :{:f}
    c    :{:f}
    i    :{:f}
    m    :{:f}
    f_t  :{:f}
    f_n  :{:f}
    wb   :{:f}
    Dl   :{:f}  
    l    :{:f}
    """.format(
        ts['y'][year],
        ts['y_s'][year],
        ts['iv'][year],
        ts['c'][year],
        ts['i'][year],
        ts['m'][year],
        ts['f_t'][year],
        ts['f_n'][year],
        ts['wb'][year],
        ts['l'][year] - ts['l'][year-1],
        ts['l'][year]
    ))
    
    # Terminate loop if system is unstable
    if ts['y'][year] < 0 or ts['y_s'][year] < 0:
        logging.critical("Unstable - exiting")
        sys.exit()

        
## %%% 3) Now distribute aggregates among micro agents and check that micro
    # and macro figures match
    temp_total_k       = 0
    temp_total_y       = 0
    temp_total_wb      = 0 
    temp_total_ft      = 0

    # Counters to tot up aggregate values (is this needed? can't we define them below as we go along?)
    agg_y        = 0 # output
    agg_iv       = 0 # inventories
    agg_k        = 0 # kapital
    agg_u        = 0 # capacity utilisation
    agg_r        = 0 # rate of profit
    agg_f_t      = 0 # gross profits
    agg_f_n      = 0 # net profits (minus interest)
    agg_f_r      = 0 # retained profits
    agg_f_d      = 0 # distributed profits
    agg_m_f      = 0 # firms' deposits
    agg_l_writedown = 0 # the writedown to loans due to bankruptcy
    agg_m_writedown = 0 # the writedown hh deposits due to bankruptcy
    agg_kshare   = 0
    
    # We use this as a check to ensure our totals add up.
    temp_num_firms     = 0

    logging.debug("*** post-expenditure firm level calculations ****")

    xts['hedge'][year] = 0
    xts['spec'][year]  = 0 
    xts['ponzi'][year] = 0
    xts['bankrupt'][year] = 0   

    for firm in firm_ts:
        # use previous-year 'firm size' and 'capital share' and
        # currently calculated stochastic element to determine
        # mark-up and distribution of goods demand for current year
        fsize  = firm['size'][year-1]   ## KK: HASN'T THIS BEEN DONE ALREADY?
        kshare = firm['k_share'][year-1]
        
        # We know supply side accomodates, so increase in capital
        # will equal investment demand
        firm['k'][year] = firm['k'][year-1] + firm['i'][year]
        agg_k += firm['k'][year]

        # Distribute aggregate demand over firms
        firm['y'][year] = ts['y'][year] * firm['y_share'][year]
        agg_y += firm['y'][year]

        # Change in inventories is difference between
        # demand and production
        firm['iv'][year] = firm['iv'][year-1] + firm['y_s'][year] - firm['y'][year]
        agg_iv += firm['iv'][year]
        
        # capacity utilisation
        firm['u'][year] = (firm['y'][year]*v)/firm['k'][year]
        #agg_u       += firm[u]*fshare

        #temp_total_wb += firm[wb_d][year]

        # Gross profits = revenues - costs
        #               = demand   - cost of supply
        firm['f_t'][year] = firm['y'][year] - firm['wb'][year]

        # Record whether firms are in a hedge, speculative, or Ponzi position
        if firm['f_t'][year] > firm['it'][year] + (0.2 * firm['l'][year]):
            xts['hedge'][year] += 1
        elif firm['f_t'][year] >= firm['it'][year]:
            xts['spec'][year] += 1
        else:
            xts['ponzi'][year] +=1
        
        # Net profits  = gross profits - net interest payments
        firm['f_n'][year] = firm['f_t'][year] - firm['it'][year]
        
        # Retained profits: only distribute dividends
        # if not making a loss
        if firm['f_n'][year] > 0:
            firm['f_r'][year] = llambda * firm['f_n'][year]
        else:
            firm['f_r'][year] = firm['f_n'][year]

        # Distributed profits
        firm['f_d'][year] = firm['f_n'][year] - firm['f_r'][year]

        # Calculate totals for profit measures
        agg_f_t += firm['f_t'][year]
        agg_f_n += firm['f_n'][year]
        agg_f_r += firm['f_r'][year]
        agg_f_d += firm['f_d'][year]
        
        # profit rate NB. need to use current year 'capital share'
        # for summary aggregation
        
        firm['r'][year] = firm['f_t'][year]/firm['k'][year]
        agg_r += firm['r'][year] * (firm['k'][year]/ ts['k'][year])
        

        # change in liquidity =
        # sources (retained earnings + CHANGE in loans) -
        # uses   (investment)
        firm['m_f'][year] = firm['m_f'][year-1] + (
            firm['f_r'][year] +
            (firm['l'][year] - firm['l'][year-1])) - (           # sources
            firm['i'][year])                                     # uses

        # BANKRUPTCY
        if firm['m_f'][year] < 0:
            # Firm is bankrupt!

            # Write off the firms loans and 'overdrafts' (negative deposits)
            # HH take the hit for now
            agg_l_writedown += firm['l'][year]
            agg_m_writedown += (firm['l'][year] - firm['m_f'][year])

            firm['m_f'][year] = 0
            firm['l'][year] = 0

            #firm[k][year] = 0
            #firm[iv][year] = 0         

            xts['bankrupt'][year] +=1

        # Can only calculate the financial variables here, after
        # distribution has taken place.
        # Deposits held by firms
        agg_m_f += firm['m_f'][year]

        # Now calculate updated capital shares
        # work out the proportion of total capital this firm has
        firm['k_share'][year] = firm['k'][year]/ts['k'][year]
        agg_kshare += firm['k_share'][year]

        # normalise this so that 1 is the mean
        # to get a 'size' measure
        firm['size'][year] = firm['k_share'][year] * num_firms 

        logging.log(8, """Firm -
        k:     {:f},
        y:     {:f},
        wb_d:  {:f},                
        share: {:f},
        size:  {:f},
        g_i:   {:f},
        i_d:   {:f},
        u:     {:f},
        r:     {:f},
        l_d:   {:f},
        m_f:   {:f},
        iv:    {:f},
        
        """.format(
            firm['k'][year],
            firm['y'][year],
            firm['wb'][year],                       
            firm['k_share'][year],
            firm['size'][year],
            firm['g_i'][year],
            firm['i'][year],
            firm['u'][year],
            firm['r'][year],
            firm['l'][year],                        
            firm['m_f'][year],
            firm['iv'][year]
        ))

    # Now back to macro to do retained profits and monetary
    # variables
    r_l = r_m = r_l_bar  ## redundant?
    ts['f_d'][year] = agg_f_d
    ts['f_r'][year] = ts['f_n'][year] - agg_f_d
    ts['m_f'][year] = agg_m_f
    ts['y_hr'][year] = ts['wb'][year] +\
        (r_m * ts['m_h'][year-1]) + agg_f_d
    ts['m_h'][year] = ts['m_h'][year-1] + \
        ts['y_hr'][year] - ts['c'][year] - agg_m_writedown # households deposits
    ts['l'][year] -= agg_l_writedown # adjust total loans for writedowns
    ts['m_t'][year] = ts['m_h'][year] + ts['m_f'][year] # total deposits
            
    logging.debug(""" Aggregated firm level numbers:
    agg_m: {}
    agg_l: {}
    agg_i: {}
    agg_y_s: {}
    agg_wb: {}
    kshare : {}
    agg_l_writedown : {}
    agg_m_writedown : {}
    """.format(
        agg_m,
        agg_l,
        agg_i,
        agg_y_s,
        agg_wb,
        agg_kshare,
        agg_l_writedown,
        agg_m_writedown,
    ))

# =============================================================================
#     logging.debug("macro k: {} micro k: {}".format(sol['k'], agg_k))
#     logging.debug("macro y: {} micro y: {}".format(sol['y'], agg_y))
#     logging.debug("macro iv:{} micro iv: {}".format(sol['iv'], agg_iv))
#     logging.debug("macro f_t: {} micro f_t: {}".format(sol['f_t'], agg_f_t))
#     logging.debug("micro f_r: {}".format(agg_f_r))
#     logging.debug("micro f_d: {}".format(agg_f_d))
#     logging.debug("micro m_f: {}".format(agg_m_f))
#     logging.debug("macro l: {}, macro m_t: {}".format(ts['l'][year], ts['m_t'][year]))
#     logging.debug("macro r: {} micro r: {}".format(sol['r'], agg_r))
#     #loging.debugg("macro u: {} micro u: {}".format(sol['u'], agg_u))
#     logging.debug("macro wb: {} micro wb: {}" .format(sol['wb'], agg_wb))
# =============================================================================
    logging.debug("*** End of year {:d} ****".format(year))

logging.info("simulation complete")

##########################################################################
######### %% END OF SIMULATION LOOP###############################
########################################################################

# Stop the timer
end_time = time.time()

# Calculate the elapsed time
elapsed_time = end_time - start_time

# Print the elapsed time in seconds
print("Elapsed time: {:.2f} seconds".format(elapsed_time))

##########################################################################
# Generate summary graphics and animations
########################################################################
# Define function that calculates variables as percent of GDP
def pct_gdp(code, ix_first, ix_last):
    return (ts[code].iloc[ix_first:ix_last] / ts['y'].iloc[ix_first:ix_last]) * 100
    

# Define function that makes charts
def mk_charts(start_year, end_year, num_years, years,
              ts, firm_ts):

    # These start at zero instead of 
    # the correct year number (why?)
    ix_first = 50         
    ix_last  = num_years  
    year     = end_year - 1
    
    
    # Money supply
    pyplot.clf()
    pyplot.plot(years[ix_first:],
                pct_gdp('m_h', ix_first, ix_last),
                label="Household deposits")
    
    pyplot.plot(years[50:],
                pct_gdp('m_f', ix_first, ix_last),
                label="Firm deposits (% GDP)")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('m.png')
#   pyplot.savefig('m.pdf')   


#    # Financial balances
#    pyplot.clf()
#    pyplot.plot(years[ix_first:],
#                map(lambda x, y: x/y,
#                    ts[s_h]._list[ix_first:ix_last],
#                    ts[y]._list[ix_first:ix_last]),
#                label="Household saving")
#    
#    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#                  ncol=2, mode="expand", borderaxespad=0.)
#    
#    pyplot.savefig('s_h.png')
#    pyplot.savefig('s_h.pdf')   
#    
    
    # Inventories
    pyplot.clf()
    pyplot.plot(years[50:],
                pct_gdp('iv', ix_first, ix_last),
                label="Inventories (% GDP)")
        
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('iv.png')
#    pyplot.savefig('iv.pdf')    
    
    
    # Investment Growth
    pyplot.clf()
    pyplot.plot(years[50:],
                ts['g_i'].iloc[ix_first:ix_last],
                label = "Investment growth")
    
#    pyplot.plot(years[50:],
#                ts[g_s]..iloc[ix_first:ix_last],
#                label = "Saving growth")
    
    pyplot.plot(years[50:],
                ts['g_y'].iloc[ix_first:ix_last],
                label = "GDP growth")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('d_i.png')
  #  pyplot.savefig('d_i.pdf')   
    
    
    # Finance categorisation
    pyplot.clf()
    pyplot.plot(years[50:],
                xts['hedge'].iloc[ix_first:ix_last],
                label = "Hedge")
    
    pyplot.plot(years[50:],
                xts['spec'].iloc[ix_first:ix_last],
                label = "Speculative")
    
    pyplot.plot(years[50:],
                xts['ponzi'].iloc[ix_first:ix_last],
                label = "Ponzi")

    pyplot.plot(years[50:],
                xts['bankrupt'].iloc[ix_first:ix_last],
                label = "Bankruptcies")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fin.png')
    #pyplot.savefig('fin.pdf')   
    
    
    # GDP Components
    pyplot.clf()
    pyplot.plot(years[50:],
                pct_gdp('i', ix_first, ix_last),
                label = "Investment % of GDP")
    
    
    pyplot.plot(years[50:],
                pct_gdp('c', ix_first, ix_last),
                label = "Consumption % of GDP")
    
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('gdp.png')
   # pyplot.savefig('gdp.pdf')   
    

    # Wage bill
    pyplot.clf()
    pyplot.plot(years[50:],
                pct_gdp('wb', ix_first, ix_last),
                label = "Wage share")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('wage_share.png')
    #pyplot.savefig('wage_share.pdf')    
    

    # Profits
    pyplot.clf()
    pyplot.plot(years[50:],
                pct_gdp('f_r', ix_first, ix_last),
                label = "Retained profits % of GDP")
    
    pyplot.plot(years[50:], 
                pct_gdp('f_d', ix_first, ix_last),
                label = "Distributed profits % of GDP")
    
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('profits.png')
    #pyplot.savefig('profits.pdf')


    # Aggregate profit share
    pyplot.clf()
    pyplot.plot(years[50:],
                pct_gdp('f_r', ix_first, ix_last),
                label = "Retained profits % of GDP")
    
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('profits.png')
    
    # Capacity utilisation
    pyplot.clf()
    pyplot.plot(years[50:], 
                ts['u'].iloc[ix_first:ix_last],
                label = "Capacity utilisation")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('u.png')
    #pyplot.savefig('u.pdf') 
    
    
    # Rate of profit
    pyplot.clf()
    pyplot.plot(years[50:], 
                ts['r'].iloc[ix_first:ix_last],
                label = "Rate of profit")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('r.png')
   # pyplot.savefig('r.pdf') 
    
    
    # Total deposits and loans
    pyplot.clf()
    pyplot.plot(years[50:], 
                ts['l'].iloc[ix_first:ix_last],
                label = "loans")
    
    pyplot.plot(years[50:], 
                ts['m_t'].iloc[ix_first:ix_last],
                label = "deposits")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('l_d.png')
   # pyplot.savefig('l_d.pdf') 


    pyplot.clf()
    pyplot.plot(years[50:],
                pct_gdp('l', ix_first, ix_last),
                label = "loans (% GDP)")
    
    pyplot.plot(years[50:],
                pct_gdp('m_t', ix_first, ix_last),
                label = "deposits (% GDP)")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('l_d.png')
    #pyplot.savefig('l_d.pdf') 
    
    
    # Firm size distribution
#    fdist_1 = []
#    fdist_2 = []
#    fdist_3 = []
#    fdist_4 = []
#    fdist_5 = []
#    fdist_6 = []
#
    
    firm_size_list   = []
    firm_k_list      = []   
    firm_r_list      = []
    firm_u_list      = []
    firm_l_d_list    = []
    firm_m_f_list    = []
    firm_lev_list    = []
    firm_m_f_k_list  = []
    firm_l_d_k_list  = []
    firm_g_list      = []
    firm_m_list      = []       
    
    for firm in firm_ts:
#        fdist_1.append(firm['size'][1])
#        fdist_2.append(firm['size'][200])
#        fdist_3.append(firm['size'][400])
#        fdist_4.append(firm['size'][600])
#        fdist_5.append(firm['size'][800])
#        fdist_6.append(firm['size'][1000])      
#    
        firm_size_list.append(firm['size'][year])
        firm_k_list.append(firm['k'][year])       
        firm_r_list.append(firm['r'][year])
        firm_u_list.append(firm['u'][year])
        firm_l_d_list.append(firm['l'][year])
        firm_m_f_list.append(firm['m_f'][year])   
        firm_m_f_k_list.append(firm['m_f'][year]/firm['k'][year])
        firm_l_d_k_list.append(firm['l'][year]/firm['k'][year])       
        firm_lev_list.append(
            (firm['l'][year]-firm['m_f'][year])/firm['k'][year]
            )
        firm_g_list.append(firm['g_i'][year])
        firm_m_list.append(firm['m'][year])
    
#     pyplot.clf()
#     
#     n, bins, patches = pyplot.hist([fdist_1,
#                                     fdist_2,
#                                     fdist_3,
#                                     fdist_4,
#                                     fdist_5,
#                                     fdist_6],
#                                    bins=100,
#                                    histtype='step')
#     
#     pyplot.savefig('fdist.png')
#     pyplot.savefig('fdist.pdf') 
# 
#     pyplot.clf()
#     
#     n, bins, patches = pyplot.hist([fdist_1,
#                                     firm_size_list],
#                                    bins=100,
#                                    histtype='step')
#     
#     pyplot.savefig('fdist-simple.png')
#     pyplot.savefig('fdist-simple.pdf')  
# 
    
    # Distribution of profit rate among firms
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_r_list,
                'r.',
                label = "Distribution of profit rates")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-r.png')
    #pyplot.savefig('fdist-r.pdf')   
    
    
    # Distribution of loans among firms
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_l_d_list,
                'r.',
                label = "Distribution of loans")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-ld.png')
    #pyplot.savefig('fdist-ld.pdf')  
    

    # Distribution of deposits among firms
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_m_f_list,
                'r.',
                label = "Distribution of deposits among firms")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-mf.png')
    #pyplot.savefig('fdist-mf.pdf') 
    
    # Distribution of profit margins among firms
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_m_list,
                'r.',
                label = "Distribution of profit margins among firms")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-m.png')
    #pyplot.savefig('fdist-mf.pdf') 
    
    
    # Distribution of net leverage among firms
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_lev_list,
                'r.',
                label = "Distribution of net leverage among firms")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-lq.png')
    #pyplot.savefig('fdist-lq.pdf')  
    

    # Detailed distribution of net leverage among firms
    # Add capital distribution and net worth??
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_l_d_k_list,
                'r.',
                label = "Loan distribution")
    pyplot.plot(firm_size_list, 
                firm_m_f_k_list,
                'g.',
                label = "Deposit distribution")
    pyplot.plot(firm_size_list, 
                firm_lev_list,
                'b.',
                label = "Net leverage among firms")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-lqd.png')
   # pyplot.savefig('fdist-lqd.pdf') 
    
    # Distribution of capacity utilisation
    # within firms
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_u_list,
                'r.',
                label = "Distribution of utilisation rates")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-u.png')
   # pyplot.savefig('fdist-u.pdf')   

    # Distribution of investment growth rates
    # within firms
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_g_list,
                'r.',
                label = "Distribution of growth rates")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-g.png')
    #pyplot.savefig('fdist-g.pdf')   

    logging.info("Finished plotting graphs")


logging.basicConfig(level=logging.INFO)

##############################################################################
# ################ %%% Execute graphs function
##############################################################################
os.chdir("C:/Users\Karsten Kohler/Dropbox/ABM project with Jo/images_pandas")
mk_charts(start_year, end_year, num_years, years,
          ts, firm_ts)


###############################################################################
###  %%% Some summary stats
###############################################################################
# Obtain firm values from last period
firm_size_last   = []
firm_k_last      = []   
firm_r_last      = []
firm_u_last      = []
firm_l_d_last    = []
firm_m_f_last    = []
firm_lev_last    = []
firm_m_f_k_last  = []
firm_l_d_k_last  = []
firm_g_last      = []
firm_m_last      = []       
    
for firm in firm_ts:
        firm_size_last.append(firm['size'][year])
        firm_k_last.append(firm['k'][year])       
        firm_r_last.append(firm['r'][year])
        firm_u_last.append(firm['u'][year])
        firm_l_d_last.append(firm['l'][year])
        firm_m_f_last.append(firm['m_f'][year])   
        firm_m_f_k_last.append(firm['m_f'][year]/firm['k'][year])
        firm_l_d_k_last.append(firm['l'][year]/firm['k'][year])       
        firm_lev_last.append(
            (firm['l'][year]-firm['m_f'][year])/firm['k'][year]
            )
        firm_g_last.append(firm['g_i'][year])
        firm_m_last.append(firm['m'][year])

# Obtain firm values from first period (after initialisation)
firm_size_first   = []
firm_k_first      = []   
firm_r_first      = []
firm_u_first      = []
firm_l_d_first    = []
firm_m_f_first    = []
firm_lev_first    = []
firm_m_f_k_first  = []
firm_l_d_k_first  = []
firm_g_first      = []
firm_m_first      = []       
    
for firm in firm_ts:
        firm_size_first.append(firm['size'][2])
        firm_k_first.append(firm['k'][2])       
        firm_r_first.append(firm['r'][2])
        firm_u_first.append(firm['u'][2])
        firm_l_d_first.append(firm['l'][2])
        firm_m_f_first.append(firm['m_f'][2])   
        firm_m_f_k_first.append(firm['m_f'][2]/firm['k'][2])
        firm_l_d_k_first.append(firm['l'][2]/firm['k'][2])       
        firm_lev_first.append(
            (firm['l'][year]-firm['m_f'][2])/firm['k'][2]
            )
        firm_g_first.append(firm['g_i'][2])
        firm_m_first.append(firm['m'][2])

# Average firm profit share 
print((sum(firm_m_last)/ num_firms)*100)

# Actual profit share (= 1 - WB/Y)
print(100 - 100*ts['wb'].iloc[year-1] / ts['y'].iloc[year-1])
print(100*ts['m'].iloc[year-1])
prof_share_agg=ts['m'] # save as series

# Calculate average and size weighted profit share per year 
# (there must be a smarter way of coding this)
profshare_avr = pd.Series(0,index=years)
profshare_wavr = pd.Series(0,index=years)
m_sum_year= pd.Series(0,index=years) 

for year in years[1:]:
    for firm in firm_ts:  
        m_sum_year[year] += firm['m'][year]
        profshare_wavr[year] += firm['m'][year]*firm['y_share'][year]

    profshare_avr[year]= m_sum_year[year]/num_firms

# size weighted profit share should be equivalent to the aggregate profit 
print(profshare_wavr[year])
print(prof_share_agg[year])
# practically yes (but there are small differences)

# =============================================================================
# ## Construct granular residual a la Gabaix 2011 (would have to do this based on top 10 firms or so)
# gran_resid = prof_share_agg - profshare_avr
# 
# ## Check explanatory power of granular residual
# import statsmodels.api as sm    
# 
# # Drop rows with missing values
# prof_share_agg = prof_share_agg.dropna()
# gran_resid =gran_resid.dropna()
# 
# # Fit OLS model
# model = sm.OLS(prof_share_agg, gran_resid)
# results = model.fit()
# Display regression results
#print(results.summary())
# =============================================================================

## Plots
import matplotlib.pyplot as plt

# Plot average firm profit share over time
plt.plot(profshare_avr[50:], linestyle='-', linewidth=2, color='k')
plt.title("Firm average profit share")
plt.show()
# is falling!

# Plot aggregate profit share over time
plt.plot(prof_share_agg[50:], linestyle='-', linewidth=2, color='k')
plt.title("Aggregate profit share")
plt.show()

# =============================================================================
# # Plot granular residual
# plt.plot(gran_resid[50:], linestyle='-', linewidth=2, color='k')
# plt.title("Granular residual")
# plt.show()
# =============================================================================

# Distribution of firm profit shares
plt.hist(firm_m_first, bins=10, edgecolor='black') 
plt.title("Initial distribution of profit shares")
plt.show()

plt.hist(firm_m_last, bins=10, edgecolor='black') 
plt.title("Final distribution of profit shares")
plt.show()

#######################################################
############## Check analytical solution
#####################################################
# Save deposit to capital ratio
dep_cap_h=ts['m_h']/ts['k']
print(dep_cap_h[year])

# Save rate of capacity utilisation
cap_util=ts['y']/ts['k']
print(cap_util[year])

# Compare different ways of calculating the profit rate
print(cap_util[year] - ts['wb'][year]/ts['k'][year])
print(ts['f_t'][year]/ts['k'][year])
print(prof_share_agg[year]*cap_util[year])
# --> all give same result

# Check aggregate consumption function
print(alpha1*((1-prof_share_agg[year])* cap_util[year] \
      + dep_cap_h[year]*r_l_bar) + alpha2*dep_cap_h[year])
print(ts['c'][year]/ts['k'][year])
# almost matches perfectly

# =============================================================================
# # Check consumption function using aggregate wage bill
# print(alpha1*((ts['wb'][year]/ts['k'][year])+ dep_cap_h[year]*r_l_bar) \
#        + alpha2*dep_cap_h[year])
# print(ts['c'][year]/ts['k'][year])
# 
# # Check consumption function using definitions
# print(alpha1*((ts['wb'][year]/ts['k'][year]) + \
#     r_m * (ts['m_h'][year]/ts['k'][year])) + alpha2*(ts['m_h'][year]/ts['k'][year]))

# =============================================================================

# Check aggregate investment function
print(gamma0 + gamma_r*((prof_share_agg[year])*cap_util[year]) + gamma_u*cap_util[year])
print(ts['i'][year]/ts['k'][year])
# matches roughly but there is a mismatch

# Check dynamic equation for cap util for t=4001
t_test=4000
print(gamma0+cap_util[t_test]*(gamma_r*prof_share_agg[t_test]+gamma_u+alpha1*(1-prof_share_agg[t_test])) \
    + dep_cap_h[t_test]*(alpha1*r_l_bar+alpha2))
print(cap_util[t_test])
# almost matches perfectly

# =============================================================================
# # Use aggregate investment instead
# print(ts['i'][4000]/ts['k'][4000] + alpha1*((1-prof_share_agg[4000])* cap_util[4000] \
#       + dep_cap_h[4000]*r_l_bar) + alpha2*dep_cap_h[4000])
# print(cap_util[4001])
# ### -> matches slightly better
# 
# =============================================================================

# Analytical solution rate of cap util
cap_util_eq=(gamma0+(alpha1*r_l_bar+alpha2)*dep_cap_h[year]) / \
(1-gamma_r*prof_share_agg[year]-gamma_u-alpha1*(1-prof_share_agg[year]))
print(cap_util_eq)

# Compare with numerical solution rate of cap util
print(cap_util[year])
# almost matches 

# Analytical solution deposit-to-capital ratio
dep_cap_h_eq= ((1-prof_share_agg[year])*cap_util[year]*(1-alpha1)) / \
(gamma0 + cap_util[year]*(gamma_r*(1-prof_share_agg[year]) + gamma_u ) \
 - r_l_bar*(1-alpha1) + alpha2)
print(dep_cap_h_eq)     
print(dep_cap_h[year])   
# almost matches  

# Construct Jacobian matrix at the equilibrium 
J11= gamma_r*prof_share_agg[year] + gamma_u + alpha1*(1-prof_share_agg[year])
J12= r_l_bar*alpha1 + alpha2
J21 = ((1-prof_share_agg[year])*(1-alpha1)*(1+gamma0) - \
       (gamma_r*prof_share_agg[year] + gamma_u)*dep_cap_h[year]*(1+r_l_bar*(1-alpha1)-alpha2)) / \
    (1+ gamma0 + gamma_r*((prof_share_agg[year])*cap_util[year]) + gamma_u*cap_util[year])**2
J22 = (1+r_l_bar*(1-alpha1)-alpha2)/(1+ gamma0 + gamma_r*((prof_share_agg[year])*cap_util[year]) + gamma_u*cap_util[year])    

J = np.array([[J11, J12], 
               [J21, J22]])
print(J)

# Obtain eigenvalues
eigenvalues, eigenvectors = np.linalg.eig(J)
print(eigenvalues)
# both real and below unity

## Confirm Jacobian by using eigenvalues to find value of cap util in period 4000
# Normalise eigenvectors
print(eigenvectors)
# Initialize an array to store the normalized eigenvectors
evecs_norm = np.copy(eigenvectors)

# Normalize the eigenvectors by dividing through by the first element
for i in range(2):
    evecs_norm[:, i] = eigenvectors[:, i] / eigenvectors[0, i]

# Print normalized eigenvectors
print(evecs_norm)

# Define the initial conditions y0
y0 = np.array([cap_util[1], dep_cap_h[1]])

# Calculate the arbitrary constants c using the normalized eigenvectors
c = np.linalg.inv(evecs_norm).dot(y0)

## Compute solution manually for cap_util at t=2000 and compare with simulated solution
t = 2000 +1
print(c[0] * eigenvalues[0] ** t +  c[1] * eigenvalues[1] ** t + cap_util_eq)
print(cap_util[t-1])
# not too different but not a perfect match either