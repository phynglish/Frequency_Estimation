def maxprob(wTpriors,datai):

    import numpy as np
    from scipy.optimize import minimize
    import pdfval

    def minfunc(wTs,datai):
        val, _, _, _, _ = pdfval.pdfval(wTs,datai)
        return -val
    
    res = minimize(minfunc, wTpriors, args = (datai), method='Nelder-Mead')
    return res.x
