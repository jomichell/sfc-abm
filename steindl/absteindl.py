import sys, time, datetime, logging
import numpy as np
import pandas as pd

from absim import SimBloc

logging.basicConfig(level=logging.INFO)

# -------------------------------------------------------------
#  start of model definition
# -------------------------------------------------------------

class Steindl(SimBloc):
    """ Main class for Steindl SFC AB model """
    
    def __init__(self, num_firms = 1000, num_periods = 500, seed = None):
        super().__init__()
        
        self.num_firms = num_firms
        self.num_periods = num_periods

        # initialise the random number generator, optionally
        # with a specific seed.
        self.rng = np.random.default_rng(seed)

    def initialise(self):
        ''' initialises firms and state variables'''
        
        logging.info("initialising model with {} periods {} firms"\
                     .format(self.num_periods, self.num_firms))

        i = self.ivars
        
        # create bank and assign state variables
        self.bank = Bank(model=self)
        self.bank[0].L   = i.L
        self.bank[0].D_h = i.D_h
        self.bank[0].D_f = i.D_h

        self[0].Y_hr = i.Y_hr
        
        # set up firms
        self.firms = Firm.init_firms(self.num_firms, model = self)

        # empty list to hold results
        self.results = []

        
    def run(self, num_periods = None):
        ''' runs the main simulation loop and collects the results '''
        
        logging.info("starting simulation for {} periods with {} firms".\
                     format(self.num_periods, self.num_firms))

        if num_periods is None:
            num_periods = self.num_periods

        for period in range(num_periods):
            logging.debug("start of period {:d}".format(period))

            # increment the lag on the state vars.
            self.incr_lag()

            # 'unpack' shortcuts to the variables
            (p, c, l1) = self.unpack()
            bank = self.bank

            # first firm-level decision bloc
            for f in self.firms:
                f.calc_production()

            # intended cons spending
            c.C_intended = p.alpha1 * l1.Y_hr + p.alpha2 * bank[-1].D_h 

            # sum investment spending and capital stock across firms
            c.update(SimBloc.aggregate(self.firms, ['I', 'K']))

            # Total intended expenditure
            c.Y_intended = c.I + c.C_intended

            # allocate random 's_share' across firms
            Firm.calc_rand_share(self.firms, self.rng)

            # second firm-level decision bloc
            for f in self.firms:
                f.calc_financing()

            # get actual spending (given possible rationing),
            # the aggregate wage bill and distributed profits
            c.update(SimBloc.aggregate(self.firms, ['Y', 'WB', 'F_d']))

            # any rationing shows up in consumption, so adjust accordingly
            c.C = c.C_intended - (c.Y_intended - c.Y)
            
            # household disposable income
            r_D = p.r_L_bar
            c.Y_hr = c.WB + (r_D * bank[-1].D_h) + c.F_d
        
            # update HH bank deposit
            bank[0].D_h += (c.Y_hr - c.C)

            # store the current variables as results
            self.update_results()

            logging.debug("end of period {:d}".format(period))
                
        logging.info("simulation complete")
        self.results = pd.DataFrame.from_dict(self.results)  

    def update_results(self):
        ''' store values of current state variables as 'results' '''

        # create dict of current state variables
        results_dict = dict(self[0])
        
        # and add state variables from bank 
        results_dict.update(dict(self.bank[0]))
        
        self.results.append(results_dict)
                
class Bank(SimBloc):
    """ Banking sector """
    
    def incr_lag(self):
        ''' The bank needs to copy last period variables as starting
        stocks for the current period '''
        super().incr_lag()
        # copy balance sheet to new period
        self.svars[0].update(self.svars[-1])
        
    def request_loan(self, req_amount):
        """ an agent has requested a loan from this bank. decide whether or how much 
        to lend """
        logging.debug("request for loan of {:.2f}".format(req_amount))
        # for now, completely elastic response:
        
        self[0].L += req_amount
        return req_amount

    def write_down(self, loan_amount, overdraft_amount):
        """ a firm has gone bust, write down the total loan and deposit
        stocks accordingly """
        logging.debug("write_down of loan {:.2f}, overdraft  {:.2f}".format(
            loan_amount, overdraft_amount))
        # A firm was 'overdrawn'. Write off the firm's loans
        # and its overdraft

        c = self[0]
        
        # write off loans
        c.L += (0 - loan_amount)

        # write off overdrafts
        c.D_f += (0 + overdraft_amount)

        # households take the hit for both
        c.D_h += (0 - (loan_amount + overdraft_amount))

class Firm(SimBloc):
    """ Firms sector """
    
    @classmethod
    def init_firms(cls, num_firms, model):
        # variables in main/macro model
        init_vars = model.ivars

        logging.debug("firm init vars from model fr: {}".format(init_vars))
        
        # create fixed number of firms with randomised shares of main variables
        # Could replace this with a dirichlet distribution?
        K_rand = model.rng.random(num_firms)
        K_rand = K_rand/K_rand.sum() 
        
        firms = []
        for K_share in K_rand:
            new_firm = Firm(model)

            new_firm.set_svars(
                lag  = 0,
                K_share = K_share, # proportion of total K
                fsize = num_firms * K_share, # firm fsize rel. to avg.
                K     = init_vars.K   * K_share,
                Y     = init_vars.Y   * K_share, #  revenue share
                u     = init_vars.u,
                r     = init_vars.r,
                D_f   = init_vars.D_f * K_share,
                L     = init_vars.L   * K_share,
                F_n   = init_vars.F_n * K_share,
                F_r   = init_vars.F_r * K_share,
                I     = init_vars.I   * K_share,
                IV    = init_vars.IV  * K_share,
            )
            firms.append(new_firm)

        return(firms)

    @classmethod
    def calc_rand_share(cls, firm_list, rng):
        """ calculates a random 'share' for each firm (replace with dirichlet?)"""
        stoch = np.array([rng.random() for f in firm_list])
        stoch = stoch/stoch.sum()

        for f, s in zip(firm_list, stoch):
            f.svars[0].s_share = s
 
    def calc_production(self):
        """ decisions made before current period revenue is known"""
        (p, c, l1) = self.unpack()
        bank = self.model.bank
        
        # Investment (expenditure) function in growth terms
        c.g_I = p.gamma0 + (p.gamma_r * l1.r) + (p.gamma_u * l1.u)

        # Investment (expenditure) in level terms
        c.I = c.g_I * l1.K
        
        # elastic supply side means increase in K will equal spending on I
        c.K = l1.K + c.I

        # Mark-up and profit margin depend on firm size
        c.tau = p.tau * l1.fsize
        c.m   = c.tau/(1+c.tau)
        
        # expected demand is previous year's demand
        # scaled by growth of capital stock
        c.Y_expected = l1.Y * (1 + c.g_I)

        # desired inventories based on expected sales
        c.IV_d = p.iota * c.Y_expected

        # Production based on expected sales and 
        # current excess inventories
        c.Y_s =  c.Y_expected + (c.IV_d - l1.IV)

        # Don't produce a negative amount
        if c.Y_s < 0:
            c.Y_s = 0
        
        # Wage bill
        c.WB   = c.Y_s * (1 - c.m)
        
        # net interest payments
        r_D = p.r_L_bar
        c.it = (p.r_L_bar * l1.L) - (r_D * l1.D_f)

        # expected current costs
        c.costs = c.WB + c.it
        
        # expected net profits
        c.F_n_e = c.Y_expected - c.costs

        # expected retained profits
        if c.F_n_e > 0:
            c.F_r_e = p.llambda * c.F_n_e
        else:
            c.F_r_e = c.F_n_e

        #  minimum desired deposits is a proportion of current costs
        c.D_f_d = p.theta * c.costs
        
        # expected deposits
        c.D_f_e = l1.D_f + c.F_r_e - c.I

        # figure out how much to borrow. if cash on hand greater
        # than desired, borrow nothing otherwise ask to borrow
        # the difference
        
        if c.D_f_e > c.D_f_d:
            L_req = 0
        else:
            L_req = c.D_f_d - c.D_f_e

        delta_L = bank.request_loan(L_req)            
        c.L = l1.L + delta_L

    def calc_financing(self):
        """ firm decisions after revenue is known plus balance sheet updates """
        (p, c, l1) = self.unpack()
        bank = self.model.bank
        m = self.model[0]

        # linear combination of capital size and stochastic 'share'
        # gives share of total intended expenditure directed to this firm
        c.Y_share = (p.zeta * l1.K_share) + ((1-p.zeta) * c.s_share)

        # expenditure intended (by consumers) for this firm
        c.Y_intended = m.Y_intended * c.Y_share

        # check if there is sufficient inventory to
        # accomodate intended expenditure

        available_stock = l1.IV + c.Y_s

        if c.Y_intended > available_stock:
            c.Y = available_stock
        else:
            c.Y = c.Y_intended
        
        # updated inventories after sales
        c.IV = available_stock - c.Y

        # consistency check for negative IV
        if c.IV < 0:
            logging.debug("IV: {}".format(c.IV))
        
        # profits
        c.F_t = c.Y - c.WB   # gross
        c.F_n = c.F_t - c.it # net

        # if profits are positive, distribute a share
        if c.F_n > 0:
            c.F_r = p.llambda * c.F_n # retained
        else:
            c.F_r = c.F_n # retained

        c.F_d = c.F_n - c.F_r # distributed

        # change in deposits = retained earnings plus change in borrowing
        # less investment spending
        delta_D_f = c.F_r + (c.L - l1.L) -  c.I

        # adjust this firm's local financial records
        c.D_f = l1.D_f + delta_D_f
            
        # adjust the bank's aggregate balance sheet
        bank[0].D_f += delta_D_f
        
        if c.D_f < 0:
            # firm is bankrupt. write off the firms loans and
            # 'overdrafts' (negative deposits)
            bank.write_down(loan_amount = c.L, overdraft_amount = 0 - c.D_f)
            c.L = 0
            c.D_f = 0

        c.r   = c.F_t / c.K # profit rate
        c.u  = c.Y * (p.v / c.K)  # utilisation. should this be Y or Y_s?
#        c.u  = c.Y_s * (p.v / c.K)  # utilisation. should this be Y or Y_s?

        # size measures for next period allocation and decision making
        c.K_share = c.K / m.K
        c.fsize = c.K_share * self.model.num_firms 

    
