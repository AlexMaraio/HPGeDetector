import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
sns.set()
import scipy.optimize as sciopt
import scipy.stats as scistats
from lmfit import Model, Parameters
import pandas as pd

def ChiSqFunc(Measured,Fitted,Errors):
    ChiSquared = 0
    for i in range(len(Measured)):
        ChiSquared += ((Measured[i] - Fitted[i])**2.0) / ((Errors[i])**2.0)
    ReducedChiSq = ChiSquared/(len(Measured)-2)
    ProbChiSq = (1.0 - scistats.chi2.cdf(ChiSquared,len(Measured)-2) )  * 100.0
    return ChiSquared, ReducedChiSq, ProbChiSq

def gauss(x,amplitude,mean,sigma,offset):
    return amplitude*np.exp(-((x-mean)**2.0)/(2.0*(sigma**2.0))) + offset

source = "Barium133-2HrRun_009_eh_1"
df = pd.read_table(source+ ".dat", sep="\s+",names = ['channel number','count number'])

df['count errors'] = np.sqrt(df['count number'])

data = df[(565<=df['channel number']) & (df['channel number']<=575)]

GaussModel = Model(gauss)
Params = Parameters()

Params.add('amplitude',value=1,vary=True)
Params.add('mean',value=570,vary=True)
Params.add('sigma',value=1,vary=True)
Params.add('offset',value=1,vary=True)

FitResult = GaussModel.fit(data['count number'],params=Params,x=data['channel number'])
print(FitResult.best_values)
FitResult.plot(yerr=data['count errors'])

thingy = ChiSqFunc(list(data['count number']),list(FitResult.best_fit),list(data['count errors']))
print(thingy)

plt.plot(data['channel number'],data['count number'])
plt.show()