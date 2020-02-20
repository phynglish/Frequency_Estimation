def buildfit(wTs,h,evecs,evals,NN):

    import numpy as np
    
    def genfuncval(tt,wTs,h,evecs,evals):
        dim = len(wTs)
        Nfreqs = int(dim/2)
        nn = (NN - 1)/2
        
        funcval = 0

        for j in range(0,dim):
            for i in range(0,Nfreqs):
                funcval = funcval + h[j]*(evecs[i,j]*np.cos(wTs[i]*(tt - nn))
                + evecs[i + Nfreqs,j]*np.sin(wTs[i]*(tt - nn)))*np.exp(-(tt -
                nn)/wTs[i + Nfreqs])/np.sqrt(evals[j])

        return funcval

    fitfunc = []
    
    for i in range(0,NN):
        funcval = genfuncval(i,wTs,h,evecs,evals)
        fitfunc.append(funcval)

    return fitfunc
