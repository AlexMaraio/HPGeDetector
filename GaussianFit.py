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

def Gauss(x,amplitude,mean,sigma,a):
    return (amplitude/(np.sqrt(2*np.pi) * sigma )) *np.exp(-((x-mean)**2.0)/(2.0*(sigma**2.0))) + a

def GaussLinear(x,amplitude,mean,sigma,a,b):
    return (amplitude/(np.sqrt(2*np.pi) * sigma ))*np.exp(-((x-mean)**2.0)/(2.0*(sigma**2.0))) + a + b*x

def GaussQuad(x,amplitude,mean,sigma,a,b,c):
    return (amplitude/(np.sqrt(2*np.pi) * sigma ))*np.exp(-((x-mean)**2.0)/(2.0*(sigma**2.0))) + a + b*x + c*x**2.0    

SetSigma = [2,3]
SetSigma = 2

BariumList = [ [[550,585,570],[1125,1140,1132],[1565,1582,1574],[1935,1960,1948],[2118,2150,2135],[2493,2526,2509],[2688,2723,2705]] , [] ]

PeakNo = 1

datadict = {'Fit type':[], 'Peak number':[], 'Peak type':[], 'Min of range':[], 'Max of range':[], 'Mean':[], 'Amplitude':[], 'Sigma':[], 'a':[], 'b':[], 'c':[], 'chisq':[], 'Reduced chisq':[], 'Probchisq':[]}

'''
We want a DataFrame object with the columns:
Fit type, peak number, peak type, min_of_range, max_of_range, mean, amplitude, sigma, a, b, c, chisq, red_chisq, probchi  
'''

for item in BariumList:
    for peak in item:
        if item == BariumList[0]:
            PeakType = "Full"
        else:
            PeakType = "Zoom"

        source = "Barium133-2HrRun_009_eh_1"
        df = pd.read_table(source+ ".dat", sep="\s+",names = ['channel number','count number'])

        df['count errors'] = np.sqrt(df['count number'])

        MinValue = peak[0]
        MaxValue = peak[1]

        data = df[(MinValue<=df['channel number']) & (df['channel number']<=MaxValue)]
        MeanValue = peak[2]
        print(MinValue,MaxValue,MeanValue)
        print(BariumList)
        
        #----------------------------------------------------------------------------------------------------------
        
        GaussModel = Model(Gauss)
        Params = Parameters()

        Params.add('amplitude',value=max(df['count number']),vary=True,min=0)
        Params.add('mean',value=MeanValue,vary=True)
        Params.add('sigma',value=1,vary=True)
        Params.add('a',value=1,vary=True)

        FitResult = GaussModel.fit(data['count number'],params=Params,x=data['channel number'])
        TempList = [ FitResult.best_values['mean'] - SetSigma *FitResult.best_values['sigma'] , FitResult.best_values['mean'] + SetSigma *FitResult.best_values['sigma'] , FitResult.best_values['mean']  ] 
        if item == BariumList[0]:
            BariumList[1].append(TempList)
            #print('MeeeeevVVV')
        print(FitResult.best_values)
        FitResult.plot(yerr=data['count errors'])

        thingy = ChiSqFunc(list(data['count number']),list(FitResult.best_fit),list(data['count errors']))
        print(thingy)

        plt.plot(data['channel number'],data['count number'])
        plt.savefig(f'Plots/Gauss/Gauss_{PeakNo}_{PeakType}')
        plt.close('all')
        #plt.show()

        datadict['Fit type'].append('Gauss')
        datadict['Peak number'].append(PeakNo)
        datadict['Peak type'].append(PeakType)
        datadict['Min of range'].append(MinValue)
        datadict['Max of range'].append(MaxValue)
        datadict['Mean'].append(FitResult.best_values['mean'])
        datadict['Amplitude'].append(FitResult.best_values['amplitude'])
        datadict['Sigma'].append(FitResult.best_values['sigma'])
        datadict['a'].append(FitResult.best_values['a'])
        datadict['b'].append(0)
        datadict['c'].append(0)
        datadict['chisq'].append(thingy[0])
        datadict['Reduced chisq'].append(thingy[1])
        datadict['Probchisq'].append(thingy[2])

        
        #----------------------------------------------------------------------------------------------------------
        
        GaussModel2 = Model(GaussLinear)
        Params2 = Parameters()

        Params2.add('amplitude',value=max(df['count number']),vary=True,min=0)
        Params2.add('mean',value=MeanValue,vary=True)
        Params2.add('sigma',value=1,vary=True)
        Params2.add('a',value=1,vary=True)
        Params2.add('b',value=1,vary=True)

        FitResult2 = GaussModel2.fit(data['count number'],params=Params2,x=data['channel number'])
        #TempList2 = [ FitResult2.best_values['mean'] - SetSigma *FitResult2.best_values['sigma'] , FitResult2.best_values['mean'] + SetSigma *FitResult2.best_values['sigma'] , FitResult2.best_values['mean']  ] 
        #BariumList[1].append(TempList2)
        print(FitResult2.best_values)
        FitResult2.plot(yerr=data['count errors'])

        thingy2 = ChiSqFunc(list(data['count number']),list(FitResult2.best_fit),list(data['count errors']))
        print(thingy2)

        plt.plot(data['channel number'],data['count number'])
        plt.savefig(f'Plots/Linear/Linear_{PeakNo}_{PeakType}')
        plt.close('all')
        #plt.show()
        
        datadict['Fit type'].append('Linear')
        datadict['Peak number'].append(PeakNo)
        datadict['Peak type'].append(PeakType)
        datadict['Min of range'].append(MinValue)
        datadict['Max of range'].append(MaxValue)
        datadict['Mean'].append(FitResult2.best_values['mean'])
        datadict['Amplitude'].append(FitResult2.best_values['amplitude'])
        datadict['Sigma'].append(FitResult2.best_values['sigma'])
        datadict['a'].append(FitResult2.best_values['a'])
        datadict['b'].append(FitResult2.best_values['b'])
        datadict['c'].append(0)
        datadict['chisq'].append(thingy2[0])
        datadict['Reduced chisq'].append(thingy2[1])
        datadict['Probchisq'].append(thingy2[2])

        #----------------------------------------------------------------------------------------------------------
        
        GaussModel3 = Model(GaussQuad)
        Params3 = Parameters()

        Params3.add('amplitude',value=max(df['count number']),vary=True,min=0)
        Params3.add('mean',value=MeanValue,vary=True)
        Params3.add('sigma',value=1,vary=True)
        Params3.add('a',value=1,vary=True)
        Params3.add('b',value=1,vary=True)
        Params3.add('c',value=1,vary=True)

        FitResult3 = GaussModel3.fit(data['count number'],params=Params3,x=data['channel number'])
        #TempList3 = [ FitResult3.best_values['mean'] - SetSigma *FitResult3.best_values['sigma'] , FitResult3.best_values['mean'] + SetSigma *FitResult3.best_values['sigma'] , FitResult3.best_values['mean']  ] 
        #BariumList[1].append(TempList3)
        print(FitResult3.best_values)
        FitResult3.plot(yerr=data['count errors'])

        thingy3 = ChiSqFunc(list(data['count number']),list(FitResult3.best_fit),list(data['count errors']))
        print(thingy3)

        plt.plot(data['channel number'],data['count number'])
        plt.savefig(f'Plots/Quad/Quad_{PeakNo}_{PeakType}')
        plt.close('all')
        #plt.show()

        datadict['Fit type'].append('Quad')
        datadict['Peak number'].append(PeakNo)
        datadict['Peak type'].append(PeakType)
        datadict['Min of range'].append(MinValue)
        datadict['Max of range'].append(MaxValue)
        datadict['Mean'].append(FitResult3.best_values['mean'])
        datadict['Amplitude'].append(FitResult3.best_values['amplitude'])
        datadict['Sigma'].append(FitResult3.best_values['sigma'])
        datadict['a'].append(FitResult3.best_values['a'])
        datadict['b'].append(FitResult3.best_values['b'])
        datadict['c'].append(FitResult3.best_values['c'])
        datadict['chisq'].append(thingy3[0])
        datadict['Reduced chisq'].append(thingy3[1])
        datadict['Probchisq'].append(thingy3[2])
        
        PeakNo +=1
    PeakNo = 1

print(datadict)

Fitdf = pd.DataFrame(datadict)
