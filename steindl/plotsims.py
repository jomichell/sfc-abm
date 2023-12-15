import matplotlib.pyplot as plt
import pandas as pd


def plot(sim, skip_periods = 50):

    results = sim.results[skip_periods:]
    
    plt.clf()
    plt.plot(results.I/results.Y)
    plt.show()
    
    
    plt.clf()
    plt.plot(results.WB/results.Y)
    plt.show()
    
    
    plt.clf()
    plt.plot(results.D_h/results.Y,
             label="Household deposits")
    
    plt.plot(results.D_f/results.Y,
             label="Firm deposits")
    
    plt.plot(results.L/results.Y,
             label="Loans")
    
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)
    plt.show()
    
    # Investment Growth
    plt.clf()
    plt.plot(results.I.pct_change(),
             label = "Investment growth")
    plt.plot(results.Y.pct_change(),
             label = "GDP growth")
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)
    plt.show()
    
    # GDP Components
    
    plt.clf()
    plt.plot(results.I/results.Y,
             label = "Investment % of GDP")
    
    plt.plot(results.C/results.Y,
             label = "Consumption % of GDP")
    
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)
    plt.show()
    
    
    # Wage share and distributed profits
    
    plt.clf()
    plt.plot(results.WB/results.Y,
             label = "Wage share")
    
    plt.plot(results.F_d/results.Y,
             label = "Distributed profits % of GDP")
    
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)
    plt.show()
    
    
    #
    # Firms distribution plots
    #
    
    firms = pd.DataFrame([f.to_dict() for f in sim.firms])
    
    figure, axis = plt.subplots(2, 2) 
    # size histogram 
    axis[0, 0].hist(firms.fsize, bins=100, histtype='step')
    axis[0, 0].set_title("distribution of size of firms")
    
    # profit scatter
    axis[0, 1].plot(firms.fsize, firms.r, 'r.')
    axis[0, 1].set_title("distribution of profit rates")
    
    # utlisation scatter
    axis[1, 0].plot(firms.fsize, firms.u, 'b.')
    axis[1, 0].set_title("distribution of utilisation rates")
    
    axis[1, 1].plot(firms.fsize, firms.g_I, 'g.')
    axis[1, 1].set_title("distribution of growth rates")
    
    plt.show()    
    
    figure, axis = plt.subplots(2, 2) 
    # size histogram 
    
    # loans scatter
    axis[0, 0].plot(firms.fsize, firms.L, 'r.')
    axis[0, 0].set_title("loans nominal")
    
    # deposits scatter
    axis[0, 1].plot(firms.fsize, firms.D_f, 'g.')
    axis[0, 1].set_title("deposits nominal")
    
    # net worth scatter
    axis[1, 0].plot(firms.fsize, firms.D_f - firms.L, 'b.')
    axis[1, 0].set_title("net worth nominal")
    
    # net leverage scatter
    axis[1, 1].plot(firms.fsize, firms.L / firms.K, 'r.',
                    label = "loans / K")
    axis[1, 1].plot(firms.fsize, firms.D_f / firms.K, 'g.',
                    label = "deposits / K")
    axis[1, 1].plot(firms.fsize, (firms.L - firms.D_f) / firms.K, 'b.',
                    label = "net leverage")
    
    axis[1, 1].legend()
    
    plt.show()    
    


