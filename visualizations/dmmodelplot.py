from array import array
import matplotlib.pyplot as plt
import numpy as np
import dmdd
import pylab as pl
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def expectvalues(element='xenon', Qmin=5, Qmax=100, binsize=1, v_lag=220., 
                 efficiency= dmdd.eff.efficiency_unit, 
                 exposure=2000, 
                 sigma_name='sigma_si', sigma_val=200,
                 mass=50.0, v_esc=540.0, 
                 v_rms=220.0, rho_x=0.3):
    """Returns the expected number of recoils for an associated energy.
    Auto pass everything. sigma val is cross section, sigma name is theory.
    """
    bins = (Qmax - Qmin)/binsize #number of bins graph has
    energy_lower = np.arange(Qmin, Qmax, binsize)
    energy_upper = np.arange(Qmin+binsize, Qmax+binsize, binsize)

    result = []

    for qmin, qmax in zip(energy_lower,energy_upper):
        n = dmdd.Nexpected(element, qmin, qmax, 
                           exposure, efficiency, 
                           sigma_name, sigma_val,  
                           mass=mass, v_esc=v_esc, 
                           v_lag=v_lag, v_rms=v_rms, 
                           rho_x = rho_x)
        result.append(n)
    result = np.array(result)
    # the problem I think is that v_lag has to be in the right place in dmdd.Nexpected, but it has to be a KWARG
    # if it comes after the other KWARGs. But when given a new value 

    return result 

def calc_xy(element='xenon', Qmin=5, Qmax=100, binsize=1, v_lag=220., 
                 efficiency= dmdd.eff.efficiency_unit, 
                 exposure=2000, 
                 sigma_name='sigma_si', sigma_val=200,
                 mass=50.0, v_esc=540.0, 
                 v_rms=220.0, rho_x=0.3): 
    """
    Calculates x and y for the graph, where x is the energy and energy bins and 
    y is the recoils per energy bin. Both x and y are arrays.
    """
    energy_lower = np.arange(Qmin, Qmax, binsize)
    energy_upper = np.arange(Qmin+binsize, Qmax+binsize, binsize)
    qmid = np.zeros(len(energy_lower))
    qmid = (energy_upper + energy_lower)/2
    for i in np.arange(len(energy_upper)): 
        qmid[i] = (energy_upper[i] + energy_lower[i])/2 

        
    data = expectvalues(element=element, Qmin=Qmin, Qmax=Qmax, exposure=exposure,
                        efficiency=efficiency, binsize=binsize,
                        v_lag=v_lag, sigma_name=sigma_name, sigma_val=sigma_val,
                        mass=mass, v_esc=v_esc, v_rms=v_rms, rho_x=rho_x) 
    data_array = np.asarray(data) #turns list into an array
    
    datalist = []
    datalist.append(qmid)
    datalist.append(data_array)
    np.savetxt('xenondata.txt', datalist)#float argument required, not np.ndarray
    
    datalist=np.loadtxt('xenondata.txt') 

    return datalist
    

def plot_xy(vlag_list, times, element='xenon', Qmin=5, Qmax=100, binsize=1, 
                 efficiency= dmdd.eff.efficiency_unit, 
                 exposure=2000, 
                 sigma_name='sigma_si', sigma_val=200,
                 mass=50.0, v_esc=540.0, 
                 v_rms=220.0, rho_x=0.3, residuals=False):
    
    """
    This function needs to take a list of velocities, calcuate the x and y for each item in the list
    (though x should be same) and plot each function on the same graph. Argument Vlag_list must be a list.
    """

    x0, y0 = calc_xy(v_lag=220, element=element, Qmin=Qmin, Qmax=Qmax, binsize=binsize, 
                 efficiency= efficiency, 
                 exposure=exposure, 
                 sigma_name=sigma_name, sigma_val=sigma_val,
                 mass=mass, v_esc=v_esc, 
                 v_rms=v_rms, rho_x=rho_x)

    for v_lag, time in zip(vlag_list, times):
        plt.figure()
        x, y = calc_xy(v_lag=v_lag, element=element, Qmin=Qmin, Qmax=Qmax, binsize=binsize, 
                 efficiency= efficiency, 
                 exposure=exposure, 
                 sigma_name=sigma_name, sigma_val=sigma_val,
                 mass=mass, v_esc=v_esc, 
                 v_rms=v_rms, rho_x=rho_x)

        
        plt.xlabel('Energy [keV]', fontsize = 13)
        if sigma_name is 'sigma_si':
            plt.title('Dark Matter Collision Spectrum for Spin Independent Model', fontsize = 15)
            plt.xlim((0,80))
        else:
            plt.title('Dark Matter Collision Spectrum for Anapole Model', fontsize = 15)
            plt.xlim((0,80))
            
        if residuals:
            plt.plot(x,((y-y0)/y0)*100, label = r'$%i\/\ \mathrm{  km/s} $' % (v_lag))
            plt.ylabel('Number of Residual Expected Recoil Events', fontsize = 13)
            plt.xlim((0,80))
            #m3 = MultipleLocator(3)
            
            #plt.axes().yaxis.set_minor_locator(m3)
        else:
            plt.plot(x,y, label = r'$%i\/\ \mathrm{  km/s} $' % (v_lag))
            plt.ylabel('Number of Expected Recoil Events', fontsize = 13)
            plt.ylim((0, 175))
            #ml = MultipleLocator(5)
            #m2 = MultipleLocator(25)
            #plt.axes().yaxis.set_major_locator(m2)
            #plt.axes().yaxis.set_minor_locator(ml)

            
        plt.legend(title = 'Relative Velocity of Dark Matter Particle')        
        fig = plt.gcf()

            
        #(:.0f)
        filename = '{}_vlag_{}_r{}.png'.format(sigma_name, time, residuals)
        
        fig.savefig(filename)
        plt.close(fig)        
        
        
    
        #the following just changes the colors and axis of graph
        ml = MultipleLocator(5)
        m2 = MultipleLocator(25)
        plt.axes().xaxis.set_minor_locator(ml)
        plt.axes().yaxis.set_major_locator(m2)
        plt.axes().yaxis.set_minor_locator(ml)

        

        #plt.show()
  
        
def calc_vlag(t, v_sun = 220., v_earth= 30.):
    """Calculates the combined velocity of the Sun and the Earth's movement around the galactic center
    as a function of days (t). t = 0 corresponds to the maximum speed v = 250. Can take a single time or multiple times.
    If given multiple times, will return as an array of velocities.
    """
    
    v_lag = 220 + 30*np.cos((2*np.pi*np.asarray(t))/(365))
    return v_lag
    
    
    

def modelplot(t1, t2, plotnumber = 10, element='xenon', Qmin=5, Qmax=100, binsize=1, v_lag=220., 
                 efficiency= dmdd.eff.efficiency_unit, 
                 exposure=2000, 
                 sigma_name='sigma_si', sigma_val=200,
                 mass=50.0, v_esc=540.0, 
                 v_rms=220.0, rho_x=0.3, residuals=False):
    
    times = np.linspace(t1, t2, plotnumber)
    #print times
    #linspace creates an array from the two given numbers 
    #including the first and the last
    #and splits them evenly based on the third number
    
    relative_v = calc_vlag(times)
    #print relative_v
    #print type(times)
    #the array from times is boolean?, so we need to make it a list to pass to plot_xy

    plot_xy(vlag_list = relative_v.tolist(), times=times, element=element, Qmin=Qmin,
                        Qmax=Qmax, exposure=exposure,
                        efficiency=efficiency, binsize=binsize,
                        sigma_name=sigma_name, sigma_val=sigma_val,
                        mass=mass, v_esc=v_esc, v_rms=v_rms, rho_x=rho_x, residuals=residuals)

def dayplot(t1, t2, element='xenon', Qmin=5., Qmax=100., binsize=10, v_lag=220., 
                 efficiency= dmdd.eff.efficiency_unit, 
                 exposure=2000, 
                 sigma_name='sigma_si', sigma_val=200,
                 mass=50.0, v_esc=540.0, 
                 v_rms=220.0, rho_x=0.3):
    plt.figure(figsize=(8,8))

    days = np.linspace(t1, t2, binsize)

    vdays = calc_vlag(days)
    
    n_expect0 = expectvalues(v_lag = 220)
    n0 = sum(n_expect0)
    print n0
    
    n_list = []
    
    for v_lag in vdays:
        n_expect = expectvalues(v_lag = v_lag)
        n_sum = sum(n_expect) 
        n_list.append(n_sum)
    plt.title('Residual Dark Matter Collisions Over Time', fontsize = 15)
    plt.xlabel('Time [Days]', fontsize = 13)
    plt.ylabel('Average Number of Residual Recoil Events', fontsize = 13)
    plt.xlim((0,400))
    plt.ylim((-40, 40))
    
    #no difference between sigma_anapole and sigma_si


    
    plt.plot(days, n_list - n0)
    plt.savefig('n_events_vs_days.png')

def duoplot(Qmin, Qmax):


    """
    This function returns both the SI and Anapole models for an energy range.
    """

    x1, y1 = calc_xy(element='xenon', Qmin=Qmin, Qmax=Qmax, binsize=1, v_lag=220., 
                 efficiency= dmdd.eff.efficiency_unit, 
                 exposure=2000, 
                 sigma_name='sigma_si', sigma_val=200,
                 mass=50.0, v_esc=540.0, 
                 v_rms=220.0, rho_x=0.3)
    x2, y2 = calc_xy(element='xenon', Qmin=Qmin, Qmax=Qmax, binsize=1, v_lag=220., 
                 efficiency= dmdd.eff.efficiency_unit, 
                 exposure=2000, 
                 sigma_name='sigma_anapole', sigma_val=200,
                 mass=50.0, v_esc=540.0, 
                 v_rms=220.0, rho_x=0.3)
    plt.figure(figsize=(8,8))
    plt.title('Title Here', fontsize = 15)
    plt.xlabel('Energy [keV]', fontsize = 13)
    plt.ylabel('Average Number of Recoil Events', fontsize = 13)
    plt.plot(x1, y1, label = 'Spin Independent', color='red')
    plt.plot(x2, y2, label = 'Anapole')
    plt.legend(title = 'Model')
    ml = MultipleLocator(5)
    m2 = MultipleLocator(10)
    plt.axes().xaxis.set_minor_locator(ml)
    plt.axes().xaxis.set_major_locator(m2)
    plt.axes().yaxis.set_minor_locator(ml)
    plt.savefig('SI_vs_Anapole.png')
                
