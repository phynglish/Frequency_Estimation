def MetHast(wTinits,wTpriors,datai,N_MH,MH_widths,Tmeas):

    import postlogprob
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    
    dim = len(wTinits)
    Nfreqs = int(dim/2)
    NN = len(datai)

    ws = np.zeros((Nfreqs,N_MH+1))
    Ts = np.zeros((Nfreqs,N_MH+1))
    
    for i in range(0,Nfreqs):
        ws[i][0] = wTinits[i]
        Ts[i][0] = wTinits[i+Nfreqs]

    initlogprob = postlogprob.postlogprob(wTinits,datai,wTpriors)
    

    accepts = 0
    
    for i in range(0,N_MH):
        wprops = []
        Tprops = []
        wcurr = []
        Tcurr = []
        for j in range(0,Nfreqs):
            wprops.append(np.random.normal(ws[j][i], MH_widths[0][j]))
            Tprops.append(np.random.normal(Ts[j][i], MH_widths[0][j+Nfreqs]))
            wcurr.append(ws[j][i])
            Tcurr.append(Ts[j][i])
        wTprops = wprops + Tprops
        prob_prop = postlogprob.postlogprob(wTprops,datai,wTpriors)
        wTcurr = wcurr + Tcurr
        prob_curr = postlogprob.postlogprob(wTcurr,datai,wTpriors)
        alpha = prob_prop - prob_curr
        a = math.log(np.random.uniform(0,1,1))
        if a < alpha:
            for j in range(0,Nfreqs):
                ws[j][i+1] = wTprops[j]
                Ts[j][i+1] = wTprops[j+Nfreqs]
            accepts = accepts + 1
        else:
            for j in range(0,Nfreqs):
                ws[j][i+1] = wTcurr[j]
                Ts[j][i+1] = wTcurr[j+Nfreqs]

    plt.figure(3)
    for i in range(0,Nfreqs):
        plt.plot((ws[i,:] - wTpriors['wmin'])/(wTpriors['wmax'] - wTpriors['wmin']))
        plt.plot((Ts[i,:] - wTpriors['Tmin'])/(wTpriors['Tmax'] - wTpriors['Tmin']))
        plt.ylim((0, 1))
        plt.show(block=False)

    Nbins = int(N_MH/100)
    fig, axs = plt.subplots(2,Nfreqs,sharey=True, tight_layout=True)
    for i in range(0,Nfreqs):
        axs[0,i].hist(ws[i,:]*(10**3)*(NN-1)/(2*math.pi*Tmeas),Nbins)
        axs[1,i].hist(Ts[i,:]*Tmeas/(NN-1),Nbins)
        plt.show(block=False)
    print('Acceptance rate was {}'.format(accepts/N_MH))
    return(ws, Ts)
