def calcSNR(hbars,datai,Nfreqs):

    import numpy as np
    
    NN = len(datai)
    msdatai = np.sum(datai*np.transpose(datai))/NN
    dim = Nfreqs*2
    m = dim + 1
    
    estsigsq = (NN/(NN - m - 2))*(msdatai - m*hbars/NN)

    SNRest = np.sqrt((m/NN)*(1 + hbars/estsigsq))

    return SNRest
