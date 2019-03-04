import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sciopt
import seaborn as sns
import pandas as pd
sns.set()
import scipy.stats as scistats


def ResolutionFit(E,a,b):
    return a/E + b

FitDF = pd.read_csv('Gamma_Peak_Stats_and_Params.csv', names=['Element', 'Fit type', 'Peak number',	'Peak type', 'Energy (keV)', 'Resolution', 'Min', 'Max', 'Mean', 'A', 'Sigma', 'Error' ,'a', 'b', 'chisq', 'Reduced chisq', 'Probchisq'],usecols=list(range(1,18)))

FitDF = FitDF[FitDF['Fit type'] == 'Linear']
FitDF = FitDF[FitDF['Peak type'] == 'Zoom']

FitDF['Energy (keV)'] = FitDF['Energy (keV)'].astype(float)
FitDF['Resolution'] = FitDF['Resolution'].astype(float).multiply(100)
SodiumDF = FitDF[(FitDF['Energy (keV)'] == 511)]
FitDF = FitDF[~(FitDF['Energy (keV)'] == 511)]


DataBa = FitDF[FitDF['Element'] == 'Barium133']
DataCo = FitDF[FitDF['Element'] == 'Cobalt60']
DataNa = FitDF[FitDF['Element'] == 'Sodium22']


Fit, Errors = sciopt.curve_fit(ResolutionFit,FitDF['Energy (keV)'],FitDF['Resolution'])
print(Fit)

TitleFont = {'size':'20', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'18'}

plt.figure(1)
Energies = np.linspace(55,1400,14000)
plt.plot(Energies,ResolutionFit(Energies,*Fit),label="Fit")

plt.errorbar(DataBa['Energy (keV)'],DataBa['Resolution'],fmt="o",color='red',label="Barium-133",markersize=7,yerr=(7.05/DataBa['Energy (keV)']),elinewidth=3,capsize=5,capthick=3,marker='o')
plt.errorbar(DataCo['Energy (keV)'],DataCo['Resolution'],fmt="o",color='blue',label="Cobalt-60",markersize=7,yerr=(7.05/DataCo['Energy (keV)']),elinewidth=3,capsize=5,capthick=3,marker='o')
plt.errorbar(DataNa['Energy (keV)'],DataNa['Resolution'],fmt="o",color='yellow',label="Sodium-22",markersize=7,yerr=(7.05/DataNa['Energy (keV)']),elinewidth=3,capsize=5,capthick=3,marker='o')
plt.errorbar(SodiumDF['Energy (keV)'],SodiumDF['Resolution'],fmt="o",color='yellow',label="",markersize=7,yerr=(7.05/SodiumDF['Energy (keV)']),elinewidth=3,capsize=5,capthick=3,marker='o')
plt.xlabel('Peak Energy [keV]',**AxTitleFont)
plt.ylabel('Energy Resolution [%]',**AxTitleFont)
plt.title('Energy Resolution as a Function of Gamma Energy',**TitleFont)
plt.legend(fontsize=16,fancybox=True,shadow=False)
#plt.show()

plt.figure(2)
plt.errorbar(DataBa['Energy (keV)'],DataBa['Resolution']-ResolutionFit(DataBa['Energy (keV)'],*Fit),fmt="o",color='red',label="Barium-133",markersize=7,yerr=(7.05/DataBa['Energy (keV)']),elinewidth=3,capsize=5,capthick=3,marker='o')
plt.errorbar(DataCo['Energy (keV)'],DataCo['Resolution']-ResolutionFit(DataCo['Energy (keV)'],*Fit),fmt="o",color='blue',label="Cobalt-60",markersize=7,yerr=(7.05/DataCo['Energy (keV)']),elinewidth=3,capsize=5,capthick=3,marker='o')
plt.errorbar(DataNa['Energy (keV)'],DataNa['Resolution']-ResolutionFit(DataNa['Energy (keV)'],*Fit),fmt="o",color='yellow',label="Sodium-22",markersize=7,yerr=(7.05/DataNa['Energy (keV)']),elinewidth=3,capsize=5,capthick=3,marker='o')
plt.errorbar(SodiumDF['Energy (keV)'],SodiumDF['Resolution']-ResolutionFit(SodiumDF['Energy (keV)'],*Fit),fmt="o",color='yellow',label="",markersize=7,yerr=(7.05/SodiumDF['Energy (keV)']),elinewidth=3,capsize=5,capthick=3,marker='o')
plt.xlabel('Peak Energy [keV]',**AxTitleFont)
plt.ylabel('Energy Resolution Residual [%]',**AxTitleFont)
plt.title('Energy Resolution Residuals',**TitleFont)
plt.legend(fontsize=16,fancybox=True,shadow=False)
plt.show()

def ChiSqFunc(Measured,Fitted,Errors,Params):
    ChiSquared = 0
    for i in range(len(Measured)):
        ChiSquared += ((Measured[i] - Fitted[i])**2.0) / ((Errors[i])**2.0)
    ReducedChiSq = ChiSquared/(len(Measured)-len(Params))
    ProbChiSq = (1.0 - scistats.chi2.cdf(ChiSquared,len(Measured)-len(Params)) )  * 100.0
    return ChiSquared, ReducedChiSq, ProbChiSq


MeasuredThungys = list(DataBa['Resolution'])+list(DataCo['Resolution']) + list(DataNa['Resolution'])
FittedThingys = list(ResolutionFit(DataBa['Energy (keV)'],*Fit)) + list(ResolutionFit(DataCo['Energy (keV)'],*Fit)) + list(ResolutionFit(DataNa['Energy (keV)'],*Fit))
ErrorThingys = list(7.05/DataBa['Energy (keV)']) + list(7.05/DataCo['Energy (keV)']) + list(7.05/DataBa['Energy (keV)'])

ChiStats = ChiSqFunc(MeasuredThungys, FittedThingys, ErrorThingys,Fit)
print(ChiStats)