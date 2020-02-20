import gendata
import pdfval
import buildpriors
import buildfit
import maxprob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time
import calcSNR
import math

##datai = []
##with open('datai36.txt', 'r') as f:
##    content = f.readlines()
##    for x in content:
##        row = x.split()
##        datai.append(float(row[0]))

coeffs = {'NN':181, 'c':0.4, 'w1':1.21, 'w2':0.226, 'w3':0.98,
          'T1':55, 'T2':65, 'T3':110, 'B1':0.3, 'B4':0.10,
          'B2':0.05, 'B5':0.22, 'B3':0.25, 'B6':0.15, 'e':0.05}

NN = coeffs['NN']
Tmeas = 60
Nfreqs = 3

wTpriors = buildpriors.buildpriors(Nfreqs,coeffs)
datai, funci = gendata.gendata(coeffs)

plt.figure(1)
plt.plot(datai)
plt.plot(funci)
plt.show(block=False)

start = time.time()
x = maxprob.maxprob(wTpriors,datai)
end = time.time()

print('Calculation time:')
print(end - start)

postval, h, hbars, evecs, evals = pdfval.pdfval(x,datai)

SNRest = calcSNR.calcSNR(hbars,datai,Nfreqs)
fitfunc = buildfit.buildfit(x,h,evecs,evals,NN)
fitfunc = fitfunc + np.mean(datai)

for i in range(0,Nfreqs):
    x[i] = x[i]*(10**3)*(NN - 1)/(2*math.pi*Tmeas)
    x[i + Nfreqs] = x[i + Nfreqs]*Tmeas/(NN - 1)

print('Max probability:')
print(postval)
print('Estimated values:')
print(x)
print('Estimated SNR:')
print(SNRest)

plt.figure(2)
plt.plot(datai)
plt.plot(fitfunc)
plt.show(block=False)
