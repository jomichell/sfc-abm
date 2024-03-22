import matplotlib.pyplot as plt
import pandas as pd


def plot(sim, skip_periods = 0):

    df = sim.results

    # skip the initial periods if required
    idx = pd.IndexSlice
    df = df.loc[idx[:, :, skip_periods:]]
    
    # ar: aggregate results; fr: firm-level results
    ar = df.xs('aggregate', level='agent').droplevel(0)

    
    # calculate some summary series
    ar['u']   = ar.Y * (sim.params.v / ar.K)
    ar['F_t'] = ar.Y - ar.WB
    ar['r']   = ar.F_t / ar.K

    # Plots
    fig, ax = plt.subplots()
    ax.plot(ar.I/ar.Y)
    ax.set_title("Investment, % of GDP")
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(ar.WB/ar.Y)
    ax.set_title("Wage share, % of GDP")
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(ar.u)
    ax.set_title("utilisation")
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(ar.r)
    ax.set_title("profit rate")
    plt.show()
    
    
    plt.plot(ar.D_h/ar.Y,
             label="Household deposits, %GDP")
    
    plt.plot(ar.D_f/ar.Y,
             label="Firm deposits, %GDP")
    
    plt.plot(ar.L/ar.Y,
             label="Loans, %GDP")
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)
    
    plt.show()
    
    # Investment Growth
    plt.clf()
    plt.plot(ar.I.pct_change(),
             label = "Investment growth")
    plt.plot(ar.Y.pct_change(),
             label = "GDP growth")
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)
    plt.show()
    
    # GDP Components
    
    plt.clf()
    plt.plot(ar.I/ar.Y,
             label = "Investment % of GDP")
    
    plt.plot(ar.C/ar.Y,
             label = "Consumption % of GDP")
    
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)
    plt.show()
    
    
    # Wage share and distributed profits
    
    plt.clf()
    plt.plot(ar.WB/ar.Y,
             label = "Wage share")
    
    plt.plot(ar.F_d/ar.Y,
             label = "Distributed profits % of GDP")
    
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)
    plt.show()



    # 
    # Firms time series plots
    #

    if sim.flags.firm_results:
        fr = df.xs('firm', level='agent')
    
        fr.I.pct_change().unstack(level='agent_idx')[1:].plot(
            title = "investment, I % growth", legend=False)
        plt.show()
        
        fr.Y.pct_change().unstack(level='agent_idx')[1:].plot(
            title = "revenue, Y % growth", legend=False)
        plt.show()
    
        fr.u.unstack(level='agent_idx')[1:].plot(
            title = "utilisation, u % ", legend=False)
        plt.show()

        fr.m.unstack(level='agent_idx')[1:].plot(
            title = "margin, m % ", legend=False)
        plt.show()

        fr.m.unstack(level='agent_idx')[1:].plot(
            title = "profit rate, r % ", legend=False)
        plt.show()
    
    #
    # Firms distribution plots
    #
    
    firms = pd.DataFrame([f.to_dict() for f in sim.firms])
    
    figure, axis = plt.subplots(2, 2) 
    # size histogram 
    axis[0, 0].hist(firms.fsize, bins=100, histtype='step')
    axis[0, 0].set_title("final distr. of size of firms")
    
    # profit scatter
    axis[0, 1].plot(firms.fsize, firms.r, 'r.')
    axis[0, 1].set_title("final distr. of r")
    
    # utlisation scatter
    axis[1, 0].plot(firms.fsize, firms.u, 'b.')
    axis[1, 0].set_title("final distr. of u")
    
    axis[1, 1].plot(firms.fsize, firms.g_I, 'g.')
    axis[1, 1].set_title("final distr. of growth")
    
    plt.show()    
    
    figure, axis = plt.subplots(2, 2) 
    # size histogram 
    
    # loans scatter
    axis[0, 0].plot(firms.fsize, firms.L, 'r.')
    axis[0, 0].set_title("final loans (nom)")
    
    # deposits scatter
    axis[0, 1].plot(firms.fsize, firms.D_f, 'g.')
    axis[0, 1].set_title("final deposits (nom)")
    
    # net worth scatter
    axis[1, 0].plot(firms.fsize, firms.D_f - firms.L, 'b.')
    axis[1, 0].set_title("net financial worth (nom)")
    
    # net leverage scatter
    axis[1, 1].plot(firms.fsize, firms.L / firms.K, 'r.',
                    label = "loans / K")
    axis[1, 1].plot(firms.fsize, firms.D_f / firms.K, 'g.',
                    label = "deposits / K")
    axis[1, 1].plot(firms.fsize, (firms.L - firms.D_f) / firms.K, 'b.',
                    label = "net leverage")
    
    axis[1, 1].legend()
    
    plt.show()    
    


