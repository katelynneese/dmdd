# y axis is number of expected events for a delta t (bin of number of days)
# x axis is time in days
# for SI and anapole models
# start at 10 day bins, that will give 36-37 points

def dayplot(t1, t2, element='xenon', Qmin=5, Qmax=100, binsize=10, v_lag=220., 
                 efficiency= dmdd.eff.efficiency_unit, 
                 exposure=2000, 
                 sigma_name='sigma_si', sigma_val=200,
                 mass=50.0, v_esc=540.0, 
                 v_rms=220.0, rho_x=0.3, residuals=False):

        days = np.linspace(t1, t2, binsize)
        #this array gives the day bins for the graph

        vday = calc_vlag(days)
        #calculates the different velocities for each day

        #plot_xy(days)
        #do I need to modify the exposure time to get what I want?
        #need y axis to be days not energy

        for v_lag in vdays:
        plt.figure()
        x = days
        y = expectvalues(element=element, Qmin=Qmin, Qmax=Qmax, exposure=exposure,
                        efficiency=efficiency, binsize=binsize,
                        v_lag=vday, sigma_name=sigma_name, sigma_val=sigma_val,
                        mass=mass, v_esc=v_esc, v_rms=v_rms, rho_x=rho_x) 
