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
        ''' sets up the model by initialising firm objects
        and state variables'''
        
        logging.info("initialising model with {} periods {} firms"\
                     .format(self.num_periods, self.num_firms))

        i = self.ivars
        
        # create HH and bank and assign state variables
        self.hh = Household(model=self)
        self.hh[0].Y_hr = i.Y_hr

        self.bank = Bank(model=self)
        self.bank[0].L   = i.L
        self.bank[0].D_h = i.D_h
        self.bank[0].D_f = i.D_h
        
        # set up firms
        self.firms = Firm.init_firms(self.num_firms, model = self)
        
    def run(self, num_periods = None):
        ''' runs the main simulation loop and collects
        the numeric results '''
        
        logging.info("starting simulation for {} periods with {} firms".\
                     format(self.num_periods, self.num_firms))

        if num_periods is None:
            num_periods = self.num_periods
            
        results = []
        
        for period in range(num_periods):
            logging.debug("start of period {:d}".format(period))

            # replaced lagged state vars with current state vars
            self.lag_svars()
            
            self.hh.lag_svars()
            self.bank.lag_svars()

            for f in self.firms:
                f.lag_svars()

            # start of decision-making
            for f in self.firms:
                f.calc_production()

            self.hh.decide_cons()
            self.calc_aggregate()
            
            # allocate random 's_share' across firms
            Firm.calc_rand_share(self.firms, self.rng)
            
            for f in self.firms:
                f.calc_financing()
 
            self.calc_aggregate2()

            #self[0].update(self.bank[0].subset(['D_f', 'D_h', 'L']))
            results.append(dict(self[0]))

            logging.debug("end of period {:d}".format(period))
            logging.debug(self.bank)

                
        logging.info("simulation complete")
        self.results = pd.DataFrame.from_dict(results)  

        
    def calc_aggregate(self):
        (p, c, l1) = self.unpack()

        # sum total investment spending across all firms
        c.update(SimBloc.aggregate(self.firms, ['I', 'K']))

        # get consumption spending from hh bloc
        c.C = self.hh[0].C
        c.Y = c.I + c.C

    def calc_aggregate2(self):
        (p, c, l1) = self.unpack()
        bank = self.bank

        # get aggregate wage bill and distributed profits
        c.update(SimBloc.aggregate(self.firms, ['WB', 'F_d']))

        r_m = p.r_l_bar

        # household disposable income
        Y_hr = c.WB + (r_m * bank[-1].D_h) + c.F_d
        
        # household disposable income
        self.hh[0].Y_hr = Y_hr

        # update HH bank deposit
        bank[0].D_h += (Y_hr - c.C)
                
class Bank(SimBloc):
    """ Banking sector """
    
    def lag_svars(self):
        ''' The bank needs to copy last period variables as starting
        stocks for the current period '''
        super().lag_svars()
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

        
class Household(SimBloc):
    """ Household sector """

    def decide_cons(self):
        (p, c, l1) = self.unpack()
        bank = self.model.bank
        
        c.C = p.alpha1 * l1.Y_hr + p.alpha2 * bank[-1].D_h # cons spending

        
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
            
    def __repr__(self):
        return """state vars:{}""".format(
            self.svars)
 
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
        c.Y_e = l1.Y * (1 + c.g_I)

        # desired inventories based on expected sales
        c.IV_d = p.iota * c.Y_e

        # Production based on expected sales and 
        # current excess inventories
        c.Y_s =  c.Y_e + (c.IV_d - l1.IV)

        # Don't produce a negative amount
        if c.Y_s < 0:
            c.Y_s = 0
        
        # Wage bill
        c.WB   = c.Y_s * (1 - c.m)
        
        # net interest payments
        r_m = p.r_l_bar
        c.it = (p.r_l_bar * l1.L) - (r_m * l1.D_f)

        # predicted current costs
        c.costs = c.WB + c.it
        
        # predicted net profits
        c.F_n_e = c.Y_e - c.costs

        # predicted retained profits
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
        # gives share of total expenditure arriving as revenue at this firm
        c.Y_share = (p.zeta * l1.K_share) + ((1-p.zeta) * c.s_share)

        # units sold
        c.Y = m.Y * c.Y_share

        # inventories 
        c.IV = l1.IV + c.Y_s - c.Y

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

    
