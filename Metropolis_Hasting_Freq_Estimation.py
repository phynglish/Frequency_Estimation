import gendata
import pdflogval
import buildinits
import buildfit
import MetHast
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time
import calcSNR
import math
import statistics

##datai = []
##with open('datai36.txt', 'r') as f:
##    content = f.readlines()
##    for x in content:
##        row = x.split()
##        datai.append(float(row[0]))

N_MH = 1000
Nfreqs = 3
coeffs = {'NN':181, 'c':0.4, 'w1':1.21, 'w2':0.226, 'w3':0.98,
          'T1':55, 'T2':65, 'T3':110, 'B1':0.3, 'B4':0.10,
          'B2':0.05, 'B5':0.22, 'B3':0.25, 'B6':0.15, 'e':0.075}
wTinits = buildinits.buildinits(Nfreqs,coeffs)
datai, funci = gendata.gendata(coeffs)
NN = len(datai)
Tmeas = 60

wTpriors = {'wmin':0,'wmax':10,'wwidth':100*2*math.pi*Tmeas/((10**3)*(NN-1)),
          'Tmin':4,'Tmax':500}

wwidths = [0.003, 0.003, 0.003, 0.01, 0.01]
Twidths = [10, 10, 10, 0.1, 0.1]
MH_widths = np.zeros((1,int(2*Nfreqs)))

for i in range(0,Nfreqs):
    MH_widths[0,i] = wwidths[i]
    MH_widths[0,i+Nfreqs] = Twidths[i]

plt.figure(1)
plt.plot(datai)
plt.plot(funci)
plt.show(block=False)

start = time.time()
ws, Ts = MetHast.MetHast(wTinits,wTpriors,datai,N_MH,MH_widths,Tmeas)
end = time.time()

print('Calculation time:')
print(end - start)

x = []
for i in range(0,Nfreqs):
    x.append(statistics.median(ws[i,:]))
for i in range(0,Nfreqs):
    x.append(statistics.median(Ts[i,:]))
postlogval, h, hbars, evecs, evals = pdflogval.pdflogval(x,datai)

SNRest = calcSNR.calcSNR(hbars,datai,Nfreqs)
fitfunc = buildfit.buildfit(x,h,evecs,evals,NN)
fitfunc = fitfunc + np.mean(datai)

for i in range(0,Nfreqs):
    x[i] = x[i]*(10**3)*(NN - 1)/(2*math.pi*Tmeas)
    x[i + Nfreqs] = x[i + Nfreqs]*Tmeas/(NN - 1)

print('Max probability:')
print(postlogval)
print('Estimated values:')
print(x)
print('Estimated SNR:')
print(SNRest)

plt.figure(2)
plt.plot(datai)
plt.plot(fitfunc)
plt.show(block=False)
