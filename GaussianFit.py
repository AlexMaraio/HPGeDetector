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

def Gauss(x,A,mean,sigma,a):
    return (A/(np.sqrt(2*np.pi) * sigma )) *np.exp(-((x-mean)**2.0)/(2.0*(sigma**2.0))) + a

def GaussLinear(x,A,mean,sigma,a,b):
    return (A/(np.sqrt(2*np.pi) * sigma ))*np.exp(-((x-mean)**2.0)/(2.0*(sigma**2.0))) + a + b*x

def GaussQuad(x,A,mean,sigma,a,b,c):
    return (A/(np.sqrt(2*np.pi) * sigma ))*np.exp(-((x-mean)**2.0)/(2.0*(sigma**2.0))) + a + b*x + c*x**2.0    

def ChToEnergy(Ch,a,b,c):
    Energy = np.sqrt( ( (Ch - c + (b**2.0 / (4*a))) )/a ) - (b / (2*a))
    return Energy

SetSigma = [2,3]
SetSigma = 2

BariumList = [ [[550,585,570,81],[1125,1140,1132,161],[1565,1582,1574,223],[1935,1960,1948,276],[2118,2150,2135,303],[2493,2526,2509,356],[2688,2723,2705,384]] , [] ]
#*SodiumList = [ [ [3557,3612,3585,511] , [8920,8963,8942,1274]] , [] ]
SodiumList = [ [ [3570,3630,3600,511] , [8966,9000,8984,1274]] , [] ]
#*CobaltList = [ [ [8206,8255,8231,1173] , [9322,9373,9347,1332]] , [] ]
CobaltList = [ [ [8250,8290,8270,1173] , [9370,9415,9392,1332]] , [] ]

SourceList = ["Barium133-24HrRun_006_eh_1","Sodium22-24HrData_002_eh_1","Cobalt60-24HrData_004_eh_1"]

PlotResolution = 300

PeakNo = int(1)

datadict = {'Element':[],'Fit type':[], 'Peak number':[], 'Peak type':[], 'Energy (keV)': [], 'Resolution':[], 'Min of range':[], 'Max of range':[], 'Mean':[], 'A':[], 'Sigma':[], 'Error on mean':[], 'a':[], 'b':[], 'c':[], 'chisq':[], 'Reduced chisq':[], 'Probchisq':[]}

# Fit parameters
#*a = 7.01164642
#*b = 8.1697501

a = 2.63977565e-06
b = 7.04767648
c = -1.91414134e-01

'''
We want a DataFrame object with the columns:
Fit type, peak number, peak type, min_of_range, max_of_range, mean, amplitude, sigma, a, b, c, chisq, red_chisq, probchi  
'''
for source in SourceList:
    if source == "Barium133-24HrRun_006_eh_1":
        ElementList = BariumList
        Element = "Barium133"
    elif source == "Sodium22-24HrData_002_eh_1":
        ElementList = SodiumList
        Element = "Sodium22"
    else:
        ElementList = CobaltList
        Element = "Cobalt60"

    for item in ElementList:
        for peak in item:
            #! Initialises the source and ranges for run

            if item == ElementList[0]:
                PeakType = "Full"
            else:
                PeakType = "Zoom"

            df = pd.read_csv(source + ".dat", sep = r"\s+", names = ['channel number','count number'])

            df['count errors'] = np.sqrt(df['count number'])
            df['count errors'] = df['count errors'].replace(0,1)

            MinValue = peak[0]
            MaxValue = peak[1]
            MeanValue = peak[2]
            PeakEnergy = peak[3]

            data = df[(MinValue<=df['channel number']) & (df['channel number']<=MaxValue)]
            #?plt.plot(data['channel number'],data['count number'])
            #?plt.show()
            #----------------------------------------------------------------------------------------------------------
            #! Gaussain + Offset model only

            GaussModel = Model(Gauss)
            Params = Parameters()

            Params.add('A',value=max(data['count number']),vary=True,min=0)
            Params.add('mean',value=MeanValue,vary=True)
            Params.add('sigma',value=1,vary=True)
            Params.add('a',value=1,vary=True)

            FitResult = GaussModel.fit(data['count number'],params=Params,x=data['channel number'])
            TempList = [ FitResult.best_values['mean'] - SetSigma *FitResult.best_values['sigma'] , FitResult.best_values['mean'] + SetSigma *FitResult.best_values['sigma'] , FitResult.best_values['mean'] , PeakEnergy ] 
            if item == ElementList[0]:
                ElementList[1].append(TempList)
                #print('MeeeeevVVV')
            #print(FitResult.best_values)
            FitResult.plot(yerr=data['count errors'],xlabel='Channel Number',ylabel='Count Number',title='Gaussian + Constant Offset')

            thingy = ChiSqFunc(list(data['count number']),list(FitResult.best_fit),list(data['count errors']))
            #print(thingy)

            #?plt.plot(data['channel number'],data['count number'])
            plt.tight_layout()
            plt.savefig(f'Plots/{Element}/Gauss/Gauss_{PeakNo}_{PeakType}.png', format='png', dpi=PlotResolution)
            plt.close('all')
            #plt.show()

            CountMax1 = FitResult.best_values['A'] / ( FitResult.best_values['sigma'] * np.sqrt(2 * np.pi ) ) + FitResult.best_values['a']
            ErrorOnMean1 = FitResult.best_values['sigma'] / np.sqrt(CountMax1)

            datadict['Element'].append(Element)
            datadict['Fit type'].append('Gauss')
            datadict['Energy (keV)'].append(PeakEnergy)
            FWHM_Energy = 2 * np.sqrt(2*np.log(2)) * ChToEnergy(FitResult.best_values['sigma'],a,b,c)
            datadict['Resolution'].append( FWHM_Energy / PeakEnergy)
            #datadict['Resolution'].append( ( (2 * np.sqrt(2*np.log(2))*FitResult.best_values['sigma'] - b ) / a ) / PeakEnergy)
            datadict['Peak number'].append(PeakNo)
            datadict['Peak type'].append(PeakType)
            datadict['Min of range'].append(MinValue)
            datadict['Max of range'].append(MaxValue)
            datadict['Mean'].append(FitResult.best_values['mean'])
            datadict['A'].append(FitResult.best_values['A'])
            datadict['Sigma'].append(FitResult.best_values['sigma'])
            #*datadict['Error on mean'].append(ErrorOnMean1)
            datadict['Error on mean'].append(1)
            datadict['a'].append(FitResult.best_values['a'])
            datadict['b'].append(0)
            datadict['c'].append(0)
            datadict['chisq'].append(thingy[0])
            datadict['Reduced chisq'].append(thingy[1])
            datadict['Probchisq'].append(thingy[2])

            
            #----------------------------------------------------------------------------------------------------------
            #! Gaussain + Linear model 

            GaussModel2 = Model(GaussLinear)
            Params2 = Parameters()

            Params2.add('A',value=max(data['count number']),vary=True,min=0)
            Params2.add('mean',value=MeanValue,vary=True)
            Params2.add('sigma',value=1,vary=True)
            Params2.add('a',value=1,vary=True)
            Params2.add('b',value=1,vary=True)

            FitResult2 = GaussModel2.fit(data['count number'],params=Params2,x=data['channel number'])
            #TempList2 = [ FitResult2.best_values['mean'] - SetSigma *FitResult2.best_values['sigma'] , FitResult2.best_values['mean'] + SetSigma *FitResult2.best_values['sigma'] , FitResult2.best_values['mean']  ] 
            #print(FitResult2.best_values)
            FitResult2.plot(yerr=data['count errors'],xlabel='Channel Number',ylabel='Count Number',title='Gaussian + Linear Polynomial')

            thingy2 = ChiSqFunc(list(data['count number']),list(FitResult2.best_fit),list(data['count errors']))
            #print(thingy2)

            #?plt.plot(data['channel number'],data['count number'])
            plt.tight_layout()
            plt.savefig(f'Plots/{Element}/Linear/Linear_{PeakNo}_{PeakType}.png', format='png', dpi=PlotResolution)
            plt.close('all')
            #plt.show()

            CountMax2 = FitResult2.best_values['A'] / ( FitResult2.best_values['sigma'] * np.sqrt(2 * np.pi ) ) + FitResult2.best_values['a'] + FitResult2.best_values['b'] * FitResult2.best_values['mean']
            ErrorOnMean2 = FitResult2.best_values['sigma'] / np.sqrt(CountMax2)
            
            datadict['Element'].append(Element)
            datadict['Fit type'].append('Linear')
            datadict['Energy (keV)'].append(PeakEnergy)
            FWHM_Energy = 2 * np.sqrt(2*np.log(2)) * ChToEnergy(FitResult2.best_values['sigma'],a,b,c)
            datadict['Resolution'].append( FWHM_Energy / PeakEnergy)
            #*datadict['Resolution'].append(  ( (2 * np.sqrt(2*np.log(2))*FitResult2.best_values['sigma'] - b) / a) / PeakEnergy)
            datadict['Peak number'].append(PeakNo)
            datadict['Peak type'].append(PeakType)
            datadict['Min of range'].append(MinValue)
            datadict['Max of range'].append(MaxValue)
            datadict['Mean'].append(FitResult2.best_values['mean'])
            datadict['A'].append(FitResult2.best_values['A'])
            datadict['Sigma'].append(FitResult2.best_values['sigma'])
            #*datadict['Error on mean'].append(ErrorOnMean2)
            datadict['Error on mean'].append(1)
            datadict['a'].append(FitResult2.best_values['a'])
            datadict['b'].append(FitResult2.best_values['b'])
            datadict['c'].append(0)
            datadict['chisq'].append(thingy2[0])
            datadict['Reduced chisq'].append(thingy2[1])
            datadict['Probchisq'].append(thingy2[2])

            #----------------------------------------------------------------------------------------------------------
            #! Gaussain + Quadratic model 

            GaussModel3 = Model(GaussQuad)
            Params3 = Parameters()

            Params3.add('A',value=max(data['count number']),vary=True,min=0)
            Params3.add('mean',value=MeanValue,vary=True)
            Params3.add('sigma',value=1,vary=True)
            Params3.add('a',value=1,vary=True)
            Params3.add('b',value=1,vary=True)
            Params3.add('c',value=1,vary=True)

            FitResult3 = GaussModel3.fit(data['count number'],params=Params3,x=data['channel number'])
            #TempList3 = [ FitResult3.best_values['mean'] - SetSigma *FitResult3.best_values['sigma'] , FitResult3.best_values['mean'] + SetSigma *FitResult3.best_values['sigma'] , FitResult3.best_values['mean']  ] 
            #print(FitResult3.best_values)
            FitResult3.plot(yerr=data['count errors'],xlabel='Channel Number',ylabel='Count Number',title='Gaussian + Quadratic Polynomial')

            thingy3 = ChiSqFunc(list(data['count number']),list(FitResult3.best_fit),list(data['count errors']))
            #print(thingy3)

            #?plt.plot(data['channel number'],data['count number'])
            plt.tight_layout()
            plt.savefig(f'Plots/{Element}/Quad/Quad_{PeakNo}_{PeakType}.png', format='png', dpi=PlotResolution)
            plt.close('all')
            #plt.show()

            CountMax3 = FitResult3.best_values['A'] / ( FitResult3.best_values['sigma'] * np.sqrt(2 * np.pi ) ) + FitResult3.best_values['a'] + FitResult3.best_values['b'] * FitResult3.best_values['mean'] + FitResult3.best_values['c'] * (FitResult3.best_values['mean'])** 2.0
            ErrorOnMean3 = FitResult3.best_values['sigma'] / np.sqrt(CountMax3)

            datadict['Element'].append(Element)
            datadict['Fit type'].append('Quad')
            datadict['Energy (keV)'].append(PeakEnergy)
            FWHM_Energy = 2 * np.sqrt(2*np.log(2)) * ChToEnergy(FitResult3.best_values['sigma'],a,b,c)
            datadict['Resolution'].append( FWHM_Energy / PeakEnergy)
            #*datadict['Resolution'].append( ( (2 * np.sqrt(2*np.log(2))*FitResult3.best_values['sigma'] - b ) / a) / PeakEnergy)
            datadict['Peak number'].append(PeakNo)
            datadict['Peak type'].append(PeakType)
            datadict['Min of range'].append(MinValue)
            datadict['Max of range'].append(MaxValue)
            datadict['Mean'].append(FitResult3.best_values['mean'])
            datadict['A'].append(FitResult3.best_values['A'])
            datadict['Sigma'].append(FitResult3.best_values['sigma'])
            #*datadict['Error on mean'].append(ErrorOnMean3)
            datadict['Error on mean'].append(1)
            datadict['a'].append(FitResult3.best_values['a'])
            datadict['b'].append(FitResult3.best_values['b'])
            datadict['c'].append(FitResult3.best_values['c'])
            datadict['chisq'].append(thingy3[0])
            datadict['Reduced chisq'].append(thingy3[1])
            datadict['Probchisq'].append(thingy3[2])
            


            Delta12 = np.abs(FitResult.best_values['mean'] - FitResult2.best_values['mean'])
            ErrorDelta12 = np.sqrt( ErrorOnMean1**2.0 + ErrorOnMean2**2.0 )

            Delta23 = np.abs(FitResult2.best_values['mean'] - FitResult3.best_values['mean'])
            ErrorDelta23 = np.sqrt( ErrorOnMean2**2.0 + ErrorOnMean3**2.0 )

            #print('Delta12 ' + str(Delta12) + ' Error on 12 ' + str(ErrorDelta12))


            #print('Delta23 ' + str(Delta23) + ' Error on 23 ' + str(ErrorDelta23))


            if PeakType == "Zoom":
                if abs(ErrorDelta23) > abs(Delta23):
                    print('best fit')
                else:
                    print('not best fit')
                    
                TitleFont = {'size':'24', 'color':'black', 'weight':'bold'} 
                AxTitleFont = {'size':'16'}
                plt.figure(figsize=(10,6))
                plt.errorbar(FitResult.best_values['mean'],1,xerr=ErrorOnMean1,elinewidth=6,capsize=10,capthick=5,marker='o',markersize=20,label="Gauss")
                plt.errorbar(FitResult2.best_values['mean'],1.25,xerr=ErrorOnMean2,elinewidth=6,capsize=10,capthick=5,marker='o',markersize=20,label="GaussLinear")
                plt.errorbar(FitResult3.best_values['mean'],1.5,xerr=ErrorOnMean3,elinewidth=6,capsize=10,capthick=5,marker='o',markersize=20,label="GaussQuad")
                plt.legend(fontsize=16,markerscale=0.33)
                plt.title(str(Element) + ' Peak Number ' + str(PeakNo),**TitleFont)
                plt.xlabel('Channel Numbers',**AxTitleFont)
                #locs, labels = plt.yticks()  
                plt.yticks([1,1.25,1.5],['Guass','GuassLinear','GaussQuad'],**AxTitleFont)
                plt.tight_layout()
                plt.savefig(f'Plots/{Element}/FitComparison_Peak{PeakNo}.png', format='png', dpi=PlotResolution)
                #plt.show()
                plt.close('all')
            
            PeakNo +=1

        PeakNo = 1

#print(datadict)

Fitdf = pd.DataFrame(datadict)
Fitdf.to_csv('Gamma_Peak_Stats_and_Params.csv')


datadict['Resolution'] = [100*x for x in datadict['Resolution']]

TitleFont = {'size':'20', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'18'}

plt.plot(datadict['Energy (keV)'],datadict['Resolution'],'bo')
plt.xlabel('Peak Energy [keV]',**AxTitleFont)
plt.ylabel('Energy Resolution [%]',**AxTitleFont)
plt.title('Energy Resolution as a Function of Gamma Energy',**TitleFont)
plt.show()
