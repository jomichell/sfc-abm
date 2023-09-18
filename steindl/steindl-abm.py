# SFC ABM steindl model from this paper:
# https://www.postkeynesian.net/working-papers/1412/
#
# The code used to generate the simulations shown in the
# paper contained some accounting errors. These errors are
# corrected in this version of the code.

import matplotlib.pyplot as pyplot
from matplotlib import animation

from series import *

import numpy as np
import sys

import time
import datetime

LOG_LEVEL = "INFO"

def log(msg, level = "INFO"):
    if not level == "DEBUG" and LOG_LEVEL == "INFO":
        print(msg)
        
#
# START INITIALISATION
#

start_year = 950
end_year   = 1450
num_years  = end_year - start_year

# Number of firms
num_firms = 1000

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
v       = 4 # capital full output ratio

llambda = 0.75 # profit retention ratio
# NB lambda is a reserved word in python.

# Exogenous
r_l_bar = 0.03
tau     = 0.25 # mark-up
pr      = 1    # labour productivity

# strength of effect of capital
# share on aggregate demand allocation
zeta    = 0.9

years = range(start_year, end_year)

# Set up dict to store time series. This will be replaced
# with a pandas dataframe 

ts = dict()
for var_name in ['c', 'k', 'l', 'f_t', 'f_r', 'f_d', 'f_n',
                 'g_i', 'g_y', 'i', 'iv', 'y_s', 'y', 'y_hr', 'y_f',
                 'm', 'm_h', 'm_f', 'm_t', 'r', 'u', 'wb', 'it',
                 'r_l', 'r_m']:
    ts[var_name] = Series(start_year, end_year)

# extra series to store summary variables which are not used
# in the model

xts = dict()
for var_name in ['hedge', 'spec', 'ponzi', 'bankrupt']:
    xts[var_name] = Series(start_year, end_year)

# create list of time-series 'balance sheets' for individual firms

firm_ts = []

for firm in range(num_firms):
    bs_series = {}
    for key in ts.keys():
        # Duplicate the aggregate balance sheet
        bs_series[key]      =  Series(start_year, end_year)

    # plus some additional summary variables
    bs_series['k_share'] = Series(start_year, end_year)
    bs_series['y_share'] = Series(start_year, end_year)
    bs_series['size']  =   Series(start_year, end_year)     
    firm_ts.append(bs_series)

# Set up initial values

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

log("""Initial variable values: 
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

print("**** Alocating initial firms' balance sheets")

# Now set up the (sparse) firm balance sheets for the initial year with
# an initial capital allocation and duplicates of the aggregates for
# rate of profit and capacity utilisation (this will result in slight faulty
# investment demand decisions by firms in the first year but won't result in
# stock-flow inconsistency).

total_k_rnd = 0

for firm in firm_ts:
    # assign random numbers to firms to allocate capital
    firm['k_rnd'] = np.random.random()
    total_k_rnd += firm['k_rnd']

# now we know the total of the random numbers, we can convert to an
# allocation

sy = start_year

for firm in firm_ts:
    # work out the proportion of total capital this firm gets
    firm['k_share'][sy] = firm['k_rnd'] / total_k_rnd

    # and how this compares to the 'average' allocation, to give
    # the relative 'size' of the firm 
    firm['size'][sy] = firm['k_share'][sy] * num_firms

    # allocate an initial capital stock to the firm
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

#
# End of initialisation
#

def solve_aggregate(year):
    """ This function solves the macro model for a given
    year. Effectively this is limited to calculating total expenditure 
    """
    print("\n\*** Solving aggregate model for year {} ***\n".format(year))

    r_l = r_l_bar
    r_m = r_l
    
    c = alpha1 * ts['y_hr'][year-1] + alpha2 * ts['m_h'][year-1] # cons spending
    i = ts['i'][year]     # investment spending
    y = c + i             # total expenditure
    
    y_s = ts['y_s'][year]   # total production (supply side)
    k = ts['k'][year-1] + i # production of capital goods (demand determined)

    iv = ts['iv'][year-1] + (y_s - y)  # chg in firms' inventories
    it = (r_l * ts['l'][year-1])-(r_m * ts['m_f'][year-1]) # net interest payments

    wb = ts['wb'][year]     # Wage bill
    f_t = y - wb            # gross profits


    f_n = f_t - it          # net profits
    l = ts['l'][year]     # firms sector loan demand

    # summary variables (not used for model solution)
    g_y = (y - ts['y'][year-1])/ts['y'][year-1]  # GDP growth rate
    g_i = (ts['i'][year]/ts['i'][year-1]) - 1    # investment growth rate 
    r = f_t / k            # aggregate profit rate
    y_f = k/v              # aggregate capacity utilisation
    u = y/y_f              # capacity utilisation
    m = (y - wb)/y         # aggregate markup
    #m_t = m_h + m_f        # total deposits


    # collect the solution values in a dict
    sol = {'c'   : c,
           'i'   : i,
           'y'   : y,
           'y_s' : y_s,
           'k'   : k,
           'iv'  : iv,
           'it'  : it,
           'wb'  : wb,
           'f_t' : f_t,
           'f_n' : f_n,
#           'f_r' : f_r,
#           'f_d' : f_d,
#           'y_hr': y_hr,
           'l'   : l,
#           'm_f' : m_f,
#           'm_h' : m_h,
#           'm_t' : m_t,
           'g_y' : g_y,
           'g_i' : g_i,
           'r'   : r,
           'y_f' : y_f,
           'u'   : y/y_f,
           'm'   : m,
           'r_l' : r_l,
           'r_m' : r_m}
    
    return sol


#
#  Main loop of iterative model solution
#

for year in years[1:]:
    # Main solution loop has three steps
    # 1) firm-level decision-making based on previous-year values.
    # 2) aggregate firm level decision and used these to solve macro model.
    # 3) allocate macro magnitudes among firms according to market shares.

    print("*** Start of year {:d} ****".format(year))

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
    

    # 1) Firm level decision making 
    print("*** Firm-level decision making ****")

    # A stochastic element is used to distribute demand across firms:
    # we assign a random 'share' magnitude to each firm which is then used
    # to distributing aggregate demand and other macro quantities

    total_rnd = 0

    for firm in firm_ts:
        # first assign random numbers to firms to allocate capital
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
        # to give the share of total expenditure that this firm will recieve
        firm['y_share'][year] = (zeta * kshare) + ((1-zeta) * s_share)

        # Investment (expenditure) function in growth terms
        firm['g_i'][year] = gamma0 + (gamma_r * firm['r'][year-1]) \
            + (gamma_u * firm['u'][year-1])

        # Investment (expenditure) in level terms
        firm['i'][year] = firm['g_i'][year] * firm['k'][year-1]

        # Sum over all investment expenditures to get aggregate figure
        agg_i += firm['i'][year]

        # We assume that larger firms can impose a greater mark-up
        f_tau  = tau * fsize

        # Thus a greater gross profit margin
        firm['m'][year]      = f_tau/(1+f_tau)
        
        # Weighted average of margins gives us aggregate (is this right?)
        agg_m += firm['y_share'][year] * firm['m'][year]        

        # Firm decides how much to produce on the basis of
        # previous year's sales, intended growth of
        # capital stock and desired inventories

        # Expected demand is previous year's demand
        # scaled by growth of capital stock
        
        y_e = (firm['y'][year-1] * (1+firm['g_i'][year]))

        # desired inventories calculated on basis
        # of expected sales
        iv_d = iota * y_e
        agg_iv_d += iv_d

        # Production based on predicted sales and 
        # current excess inventories
        firm['y_s'][year] = iv_d - firm['iv'][year-1] + y_e
        agg_y_s += firm['y_s'][year]

        # Wage bill = value of production y_s less proft margin
        # times that value of production
        firm['wb'][year]   = firm['y_s'][year] - \
            (firm['m'][year] * firm['y_s'][year])
        agg_wb          += firm['wb'][year]

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

        # Predicted retained profits
        if f_n_e > 0:
            f_r_e = llambda * f_n_e
        else:
            f_r_e = f_n_e

        # Assume desired minimum liquidity is a proportion of current costs
        lq_d = theta * f_lb
        agg_lq_d += lq_d

        # Calculate expected liquidity this period on the basis of
        # expected profits

        lq_e = firm['m_f'][year-1] + f_r_e - firm['i'][year]        
        agg_lq_e += lq_e

        # Compare expected to desired liquidity and
        # decide on loan demand accordingly

        if lq_e > lq_d:
            firm['l'][year] = firm['l'][year-1]
        else:
            firm['l'][year] = firm['l'][year-1] + lq_d - lq_e

        # Sum loans outstanding to get total
        agg_l  += firm['l'][year]
        agg_Dl_d += firm['l'][year] - firm['l'][year-1]

        log("""Firm -
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
        ), "DEBUG")

    log("""Inserting aggregates from micro:
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
    
    ts['i'][year]   = agg_i
    ts['y_s'][year] = agg_y_s
    ts['wb'][year]  = agg_wb
    ts['m'][year]   = agg_m
    ts['l'][year]   = agg_l


    print(""" summary firm level aggregates:
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

    # Solve the macro model to get aggregate numbers.

    sol = solve_aggregate(year)
    
    print("""macro solutions
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
        sol['y'],
        sol['y_s'],
        sol['iv'],
        sol['c'],
        sol['i'],
        sol['m'],                
        sol['f_t'],
        sol['f_n'],      
        sol['wb'],
        sol['l'] - ts['l'][year-1],      
        sol['l']               
        ))

    if sol['y'] < 0 or sol['y_s'] < 0:
        print("Unstable - exiting")
        sys.exit()

    # 
    for key in sol.keys():
        ts[key][year] = sol[key]
        
    #
    # Now distribute aggregates among micro agents and check that micro
    # and macro figures match
    #

    temp_total_k       = 0


    temp_total_y       = 0
    temp_total_wb      = 0 
    temp_total_ft      = 0

    # Counters to tot up aggregate values

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

    agg_writedown = 0 # the writedown to loans and hh deposits due to bankruptcy
    
    agg_kshare   = 0
    
    # We use this as a check to ensure our totals add up.
    temp_num_firms     = 0

    #
    # Post-aggregate solution
    #

    print("*** post-expenditure firm level calculations ****")


    xts['hedge'][year] = 0
    xts['spec'][year]  = 0 
    xts['ponzi'][year] = 0
    xts['bankrupt'][year] = 0   


    for firm in firm_ts:
        # use previous-year 'firm size' and 'capital share' and
        # currently calculated stochastic element to determine
        # mark-up and distribution of goods demand for current year
        fsize  = firm['size'][year-1]
        kshare = firm['k_share'][year-1]
        
        # We know supply side accomodates, so increase in capital
        # will equal investment demand
        firm['k'][year]   = firm['k'][year-1] + firm['i'][year]

        agg_k          += firm['k'][year]

        # ***************
        # Distribute aggregate demand between firms
        # ***************

        # share of demand affected by relative size
        # and stochastic element. Parameter zeta determines
        # relative weights
        firm['y'][year]   = ts['y'][year] * firm['y_share'][year]
        agg_y           += firm['y'][year]

        # Change in inventories is difference between
        # demand and production
        firm['iv'][year]  = firm['iv'][year-1] + firm['y_s'][year] - firm['y'][year]
        agg_iv          += firm['iv'][year]
        
        # capacity utilisation
        firm['u'][year]   = (firm['y'][year]*v)/firm['k'][year]
        #agg_u       += firm[u]*fshare

        
        #temp_total_wb += firm[wb_d][year]

        
        # Gross profits = revenues - costs
        #               = demand   - cost of supply
        firm['f_t'][year] = firm['y'][year] - firm['wb'][year]


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
        agg_f_t        += firm['f_t'][year]
        agg_f_n        += firm['f_n'][year]
        agg_f_r        += firm['f_r'][year]
        agg_f_d        += firm['f_d'][year]
        
        # profit rate NB. need to use current year 'capital share'
        # for summary aggregation
        
        firm['r'][year]   = firm['f_t'][year]/firm['k'][year]
        agg_r          += firm['r'][year] * (firm['k'][year]/ ts['k'][year])
        

        # Can now calculate the change in liquidity =
        # sources (retained earnings + CHANGE in loans) -
        # uses    (investment)

        
        firm['m_f'][year] = firm['m_f'][year-1] + (
            firm['f_r'][year] +
            (firm['l'][year] - firm['l'][year-1])) - (           # sources
            firm['i'][year])                                     # uses

        #
        # BANKRUPTCY
        #

        #
        # Check if this is correct
        #
        
        if firm['m_f'][year] < 0:
            # Firm is bankrupt!
            firm['m_f'][year] = 0


            
            # Write off its loans
            # HH take the hit for now
            agg_writedown += (firm['l'][year] - firm['m_f'][year])
            
            #ts['m_h'][year] -= (firm['l'][year] - firm['m_f'][year])
            #ts['m_t'][year] -= (firm['l'][year] - firm['m_f'][year])
            #ts['l'][year] -= firm['l'][year]                        
            firm['l'][year] = 0
            #firm[k][year] = 0
            #firm[iv][year] = 0         

            xts['bankrupt'][year] +=1

        #
        # Can only calculate the financial variables here, after
        # distribution has taken place.
        #
            
        agg_m_f        += firm['m_f'][year]

        # Now calculate updated capital shares
        #
        # work out the proportion of total capital this firm has

        firm['k_share'][year] = firm['k'][year]/ts['k'][year]
        agg_kshare          += firm['k_share'][year]

        
        # normalise this so that 1 is the mean
        # to get a 'size' measure
        firm['size'][year] = firm['k_share'][year] * num_firms 


        log("""Firm -
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
        ), "DEBUG")

    # Now back to macro to do retained profits and monetary
    # variables
    
    r_l = r_m = r_l_bar
    
    ts['f_d'][year] = agg_f_d
    ts['f_r'][year] = ts['f_n'][year] - agg_f_d
    
    ts['m_f'][year] = agg_m_f

    ts['y_hr'][year] = ts['wb'][year] +\
        (r_m * ts['m_h'][year-1]) + agg_f_d
    
    ts['m_h'][year] = ts['m_h'][year-1] + \
        ts['y_hr'][year] - ts['c'][year] - agg_writedown # households deposits

    ts['l'][year] -= agg_writedown # adjust total loans for writedowns
    
    ts['m_t'][year] = ts['m_h'][year] + ts['m_f'][year] # total deposits

            
    log(""" Aggregated firm level numbers:
    agg_m: {}
    agg_i: {}
    agg_y_s: {}
    agg_wb: {}
    agg_l: {}
    kshare : {}
    """.format(
        agg_m,
        agg_i,
        agg_y_s,
        agg_wb,
        agg_l,
        agg_kshare
    ))

    log("macro k: {} micro k: {}".format(sol['k'], agg_k))
    log("macro y: {} micro y: {}".format(sol['y'], agg_y))
    log("macro iv:{} micro iv: {}".format(sol['iv'], agg_iv))
    log("macro f_t: {} micro f_t: {}".format(sol['f_t'], agg_f_t))
    log("micro f_r: {}".format(agg_f_r))
    log("micro f_d: {}".format(agg_f_d))
    log("micro m_f: {}".format(agg_m_f))
    log("macro l: {}, macro m_t: {}".format(sol['l'], ts['m_t'][year]))
    log("macro r: {} micro r: {}".format(sol['r'], agg_r))
    #log("macro u: {} micro u: {}".format(sol['u'], agg_u))
    log("macro wb: {} micro wb: {}" .format(sol['wb'], agg_wb))

    print("*** End of year {:d} ****".format(year))


#
# Generate summary graphics and animations
#
    
def pct_gdp(code, ix_first, ix_last):
    return [x/y for x,y in zip(ts[code]._list[ix_first:ix_last],
                               ts['y']._list[ix_first:ix_last])]

def mk_charts(start_year, end_year, num_years, years,
              ts, firm_ts):
    
    #
    # CHARTS
    #
    
    # These start at zero instead of 
    # the correct year number (why?)
    
    ix_first = 50         
    ix_last  = num_years  
    year     = end_year - 1
    
    #
    # Money supply
    #


    
    pyplot.clf()
    pyplot.plot(years[ix_first:],
                pct_gdp('m_h', ix_first, ix_last),
                label="Household deposits")
    
    pyplot.plot(years[50:],
                pct_gdp('m_f', ix_first, ix_last),
                label="Firm deposits")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('m.png')
#    pyplot.savefig('m.pdf')   


    
#    #
#    # Financial balances
#    #
#    
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
    
    #
    # Inventories
    #
    
    pyplot.clf()
    pyplot.plot(years[50:],
                pct_gdp('iv', ix_first, ix_last),
                label="Inventories")
        
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('iv.png')
#    pyplot.savefig('iv.pdf')    
    
    
    #
    # Investment Growth
    #
    
    pyplot.clf()
    pyplot.plot(years[50:],
                ts['g_i']._list[ix_first:ix_last],
                label = "Investment growth")
    
#    pyplot.plot(years[50:],
#                ts[g_s]._list[ix_first:ix_last],
#                label = "Saving growth")
    
    pyplot.plot(years[50:],
                ts['g_y']._list[ix_first:ix_last],
                label = "GDP growth")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('d_i.png')
    pyplot.savefig('d_i.pdf')   
    
    
    #
    # Finance categorisation
    #
    
    pyplot.clf()
    pyplot.plot(years[50:],
                xts['hedge']._list[ix_first:ix_last],
                label = "Hedge")
    
    pyplot.plot(years[50:],
                xts['spec']._list[ix_first:ix_last],
                label = "Speculative")
    
    pyplot.plot(years[50:],
                xts['ponzi']._list[ix_first:ix_last],
                label = "Ponzi")

    pyplot.plot(years[50:],
                xts['bankrupt']._list[ix_first:ix_last],
                label = "Bankruptcies")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fin.png')
    pyplot.savefig('fin.pdf')   
    
    
    #
    # GDP Components
    #
    
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
    pyplot.savefig('gdp.pdf')   
    
    
    
    #
    # Wage share
    #
    
    pyplot.clf()
    pyplot.plot(years[50:],
                pct_gdp('wb', ix_first, ix_last),
                label = "Wage share")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('wage_share.png')
    pyplot.savefig('wage_share.pdf')    
    
    
    #
    # Profits
    #
    
    pyplot.clf()
    pyplot.plot(years[50:], 
                ts['m']._list[ix_first:ix_last],
                label = "Gross profit margin % of GDP")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('margin.png')
    pyplot.savefig('margin.pdf')    
    
    
    #
    # Profits
    #
    
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
    pyplot.savefig('profits.pdf')

    
    
    
    #
    # Capacity utilisation
    #
    
    pyplot.clf()
    pyplot.plot(years[50:], 
                ts['u']._list[ix_first:ix_last],
                label = "Capacity utilisation")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('u.png')
    pyplot.savefig('u.pdf') 
    
    #
    # Rate of profit
    #
    
    pyplot.clf()
    pyplot.plot(years[50:], 
                ts['r']._list[ix_first:ix_last],
                label = "Rate of profit")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('r.png')
    pyplot.savefig('r.pdf') 
    
    
    #
    # Total deposits and loans
    #

    
    pyplot.clf()
    pyplot.plot(years[50:], 
                ts['l']._list[ix_first:ix_last],
                label = "loans")
    
    pyplot.plot(years[50:], 
                ts['m_t']._list[ix_first:ix_last],
                label = "deposits")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('l_d.png')
    pyplot.savefig('l_d.pdf') 


    
    pyplot.clf()
    pyplot.plot(years[50:],
                pct_gdp('l', ix_first, ix_last),
                label = "loans")
    
    pyplot.plot(years[50:],
                pct_gdp('m_t', ix_first, ix_last),
                label = "deposits")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('l_d.png')
    pyplot.savefig('l_d.pdf') 
    
    
    #
    # Firm size distribution
    #
    
    
    fdist_1 = []
    fdist_2 = []
    fdist_3 = []
    fdist_4 = []
    fdist_5 = []
    fdist_6 = []
    fdist_7 = []
    
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
    
    for firm in firm_ts:
        fdist_1.append(firm['size'][950])
        fdist_2.append(firm['size'][1000])
        fdist_3.append(firm['size'][1050])
        fdist_4.append(firm['size'][1100])
        fdist_5.append(firm['size'][1150])
        fdist_6.append(firm['size'][1200])
        fdist_7.append(firm['size'][1249])      
    
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
    
    pyplot.clf()
    
    n, bins, patches = pyplot.hist([fdist_1,
                                    fdist_2,
                                    fdist_3,
                                    fdist_4,
                                    fdist_5,
                                    fdist_6,
                                    fdist_7],
                                   bins=100,
                                   histtype='step')
    
    pyplot.savefig('fdist.png')
    pyplot.savefig('fdist.pdf') 

    pyplot.clf()
    
    n, bins, patches = pyplot.hist([fdist_1,
                                    firm_size_list],
                                   bins=100,
                                   histtype='step')
    
    pyplot.savefig('fdist-simple.png')
    pyplot.savefig('fdist-simple.pdf')  

    
    #
    # Distribution of profit rate among firms
    #


    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_r_list,
                'r.',
                label = "Distribution of profit rates")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-r.png')
    pyplot.savefig('fdist-r.pdf')   
    
    
    
    #
    # Distribution of loans among firms
    #
    
    
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_l_d_list,
                'r.',
                label = "Distribution of loans")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-ld.png')
    pyplot.savefig('fdist-ld.pdf')  
    
    
    #
    # Distribution of deposits among firms
    #
    
    
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_m_f_list,
                'r.',
                label = "Distribution of deposits among firms")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-mf.png')
    pyplot.savefig('fdist-mf.pdf')  
    

    
    #
    # Distribution of net leverage among firms
    #
    
    
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_lev_list,
                'r.',
                label = "Distribution of net leverage among firms")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-lq.png')
    pyplot.savefig('fdist-lq.pdf')  
    

    
    #
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
    pyplot.savefig('fdist-lqd.pdf') 
    


    
    #
    # Distribution of capacity utilisation
    # within firms
    #
    
    
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_u_list,
                'r.',
                label = "Distribution of utilisation rates")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-u.png')
    pyplot.savefig('fdist-u.pdf')   


    
    #
    # Distribution of investment growth rates
    # within firms
    #
    
    
    pyplot.clf()
    pyplot.plot(firm_size_list, 
                firm_g_list,
                'r.',
                label = "Distribution of growth rates")
    
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
    
    pyplot.savefig('fdist-g.png')
    pyplot.savefig('fdist-g.pdf')   

    print("**** Finished plotting graphs")

    
mk_charts(start_year, end_year, num_years, years,
          ts, firm_ts)


