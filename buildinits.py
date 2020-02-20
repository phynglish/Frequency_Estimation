def buildinits(Nfreqs,coeffs):

    ws = [coeffs['w1'], coeffs['w2'], coeffs['w3']]
    Ts = [coeffs['T1'], coeffs['T2'], coeffs['T3']]
    
    if Nfreqs <= 3:
        winits = ws[0:Nfreqs]
        Tinits = Ts[0:Nfreqs]
        wTinits = winits + Tinits
    else:
        for i in range(0,Nfreqs - 3):
            ws.append(0.5)
            Ts.append(50)
        wTinits = ws + Ts
    
    return wTinits
