def pdflogval(wTs,datai):

    import numpy as np
    import scipy as sp
    import math

    def gmatrix(wTs,NN):
        def sjksums(w1,w2,T1,T2,NN):
            nn = (NN - 1)/2
            sjk = 0
            for ii in range(0,NN):
                sjk = (sjk +
                       np.sin(w1*(ii - nn))*np.sin(w2*(ii - nn))*np.exp(-(ii - nn)*(1/T1 + 1/T2)))
            return sjk

        def cjksums(w1,w2,T1,T2,NN):
            nn = (NN - 1)/2
            cjk = 0
            for ii in range(0,NN):
                cjk = (cjk +
                       np.cos(w1*(ii - nn))*np.cos(w2*(ii - nn))*np.exp(-(ii - nn)*(1/T1 + 1/T2)))
            return cjk

        def scjksums(w1,w2,T1,T2,NN):
            nn = (NN - 1)/2
            scjk = 0
            for ii in range(0,NN):
                scjk = (scjk +
                        np.cos(w1*(ii - nn))*np.sin(w2*(ii - nn))*np.exp(-(ii - nn)*(1/T1 + 1/T2)))
            return scjk

        def cjkconstsums(w1,T1,NN):
            nn = (NN - 1)/2
            cjkconst = 0
            for ii in range(0,NN):
                cjkconst = (cjkconst + np.cos(w1*(ii - nn))*np.exp(-(ii - nn)/T1))
            return cjkconst

        def sjkconstsums(w1,T1,NN):
            nn = (NN - 1)/2
            sjkconst = 0
            for ii in range(0,NN):
                sjkconst = (sjkconst + np.sin(w1*(ii - nn))*np.exp(-(ii - nn)/T1))
            return sjkconst
        
        dim = len(wTs)
        Nfreqs = int(dim/2)
        gjk = np.zeros((dim + 1,dim + 1))

        gjk[dim,dim] = NN
        
        for j in range(0, Nfreqs):
            for i in range(0, Nfreqs):
                gjk[i,j] = cjksums(wTs[i],wTs[j],wTs[i + Nfreqs],wTs[j + Nfreqs],NN)
                gjk[i + Nfreqs,j + Nfreqs] = sjksums(wTs[i],wTs[j],wTs[i + Nfreqs],wTs[j + Nfreqs],NN)
                gjk[i + Nfreqs,j] = scjksums(wTs[j],wTs[i],wTs[j + Nfreqs],wTs[i + Nfreqs],NN)
                gjk[i,j + Nfreqs] = scjksums(wTs[i],wTs[j],wTs[i + Nfreqs],wTs[j + Nfreqs],NN)
                gjk[dim,j] = cjkconstsums(wTs[j],wTs[j + Nfreqs],NN)
                gjk[dim,j + Nfreqs] = sjkconstsums(wTs[j],wTs[j + Nfreqs],NN)
                gjk[i,dim] = cjkconstsums(wTs[i],wTs[i + Nfreqs],NN)
                gjk[i + Nfreqs, dim] = sjkconstsums(wTs[i],wTs[i + Nfreqs],NN)

        return gjk
        

    def calchbars(wTs,evals,evecs,NN):

        dim = len(wTs)
        Nfreqs = int(dim/2)
        m = dim + 1
        
        h = np.zeros(m)
        hs = np.zeros(m)
        nn = (NN - 1)/2

        for ll in range(0,NN):
            for j in range(0,m):
                for i in range(0,Nfreqs):
                    h[j] = h[j] + datai[ll]*((evecs[i,j]*np.cos(wTs[i]*(ll - nn)) +
                            evecs[i + Nfreqs,j]*np.sin(wTs[i]*(ll - nn)))*np.exp(-(ll - nn)/wTs[i + Nfreqs]) + (evecs[dim,j]/Nfreqs))/np.sqrt(evals[j])
        
        for i in range(0,m):
            hs[i] = h[i]**2

        hbars = 0
        for i in range(0,m):
            hbars = hbars + hs[i]/m

        return (hbars,h)

    NN = len(datai)
    msdatai = np.sum(datai*np.transpose(datai))/NN

    gjk = gmatrix(wTs,NN)

    evals, evecs = np.linalg.eig(gjk)

    dim = len(wTs)
    m = dim + 1
    
    pref = 1
    for i in range(0,m):
            pref = pref/np.sqrt(evals[i])
    
    hbars,h = calchbars(wTs,evals,evecs,NN)
    
    samplogdval = (NN - m)*math.log(1/np.sqrt(2*math.pi)) + math.log(pref) + (NN - 2 - m)*math.log(np.sqrt(2)) + \
                  math.log(sp.special.gamma((NN - m)/2)) + ((m - NN)/2)*math.log(NN*msdatai - m*hbars)

    return (samplogdval,h,hbars,evecs,evals)
