def gendata(coeffs):

    import numpy as np

    def genfuncval(tt,coeffs):
        funcval = (coeffs['c'] + (coeffs['B1']*np.cos(coeffs['w1']*tt)
            + coeffs['B4']*np.sin(coeffs['w1']*tt))*np.exp(-tt/coeffs['T1'])
            + (coeffs['B2']*np.cos(coeffs['w2']*tt)
            + coeffs['B5']*np.sin(coeffs['w2']*tt))*np.exp(-tt/coeffs['T2'])
            + (coeffs['B3']*np.cos(coeffs['w3']*tt)
            + coeffs['B6']*np.sin(coeffs['w3']*tt))*np.exp(-tt/coeffs['T3']))
        return funcval

    def genfunc(coeffs):
        funci = []
        for tt in range(0,coeffs['NN']):
            fv = genfuncval(tt,coeffs)
            funci.append(fv)
        
        funci = funci - np.amin(funci)
        funci = funci/np.amax(funci)

        return funci

    funci = genfunc(coeffs)

    datai = []
   
    for tt in range(0,coeffs['NN']):
        dv = funci[tt] + coeffs['e']*np.random.normal()
        datai.append(dv)
    return (datai, funci)
