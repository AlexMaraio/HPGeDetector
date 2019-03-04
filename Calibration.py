import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sciopt
import seaborn as sns
import pandas as pd
sns.set()
import scipy.stats as scistats

source = "MysterySource-2HrRun_001_eh_1"
df = pd.read_csv(source+ ".dat", sep=r"\s+",names = ['channel number','count number'])

TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'22'}


def Linear(E,a,b):
    return E*a + b

def Quadratic(E,a,b,c):
    return a*E**2.0 + E*b + c

def ChiSqFunc(Measured,Fitted,Errors,Params):
    ChiSquared = 0
    for i in range(len(Measured)):
        ChiSquared += ((Measured[i] - Fitted[i])**2.0) / ((Errors[i])**2.0)
    ReducedChiSq = ChiSquared/(len(Measured)-len(Params))
    ProbChiSq = (1.0 - scistats.chi2.cdf(ChiSquared,len(Measured)-len(Params)) )  * 100.0
    return ChiSquared, ReducedChiSq, ProbChiSq


FitDF = pd.read_csv('Gamma_Peak_Stats_and_Params.csv', names=['Element', 'Fit type', 'Peak number',	'Peak type', 'Energy (keV)', 'Resolution', 'Min', 'Max', 'Mean', 'A', 'Sigma', 'Error' ,'a', 'b', 'chisq', 'Reduced chisq', 'Probchisq'],usecols=list(range(1,18)))

FitDF = FitDF[FitDF['Fit type'] == 'Linear']
FitDF = FitDF[FitDF['Peak type'] == 'Zoom']

FitDF['Mean'] = FitDF['Mean'].astype(float)
FitDF['Error'] = FitDF['Error'].astype(float)

FitDF['Energy (keV)'] = FitDF['Energy (keV)'].astype(float)

plt.plot(FitDF['Energy (keV)'],FitDF['Mean'],'bo')

plt.show()

DataBa = FitDF[FitDF['Element'] == 'Barium133']
DataCo = FitDF[FitDF['Element'] == 'Cobalt60']
DataNa = FitDF[FitDF['Element'] == 'Sodium22']



ParamsLinear, ErrorsLinear = sciopt.curve_fit(Quadratic,FitDF['Energy (keV)'],FitDF['Mean'])
print(ParamsLinear)

rough_EC_a = a = ParamsLinear[0]
rough_EC_b = b = ParamsLinear[1]
rough_EC_c = c = ParamsLinear[2]

Energies = np.sqrt( (FitDF['Mean'] - c + (b**2.0/ (4 * a)))/a ) - (b/(2*a))
#Energies = (FitDF['Mean'] - b)/(a) # keV

plt.errorbar(DataBa['Energy (keV)'],DataBa['Mean'], fmt='ro',label="Barium-133", markersize=10.5,yerr=DataBa['Error'].values)
plt.errorbar(DataCo['Energy (keV)'],DataCo['Mean'],fmt='bo',label="Cobalt-60", markersize=10.5,yerr=DataCo['Error'].values)
plt.errorbar(DataNa['Energy (keV)'],DataNa['Mean'],fmt='yo',label="Sodium-22", markersize=10.5,yerr=DataNa['Error'].values)

plt.plot(Energies,Quadratic(Energies,*ParamsLinear),label='Fit')
plt.title('Energy calibration plot with fit',**TitleFont)
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Channel Number',**AxTitleFont)
plt.legend()
plt.show()

Residuals = []
'''
for i in range(len(ChannelNos)):
    Residuals.append(ChannelNos.iloc[i] - linear(Energies.iloc[i],*Params))
'''
plt.errorbar(DataBa['Energy (keV)'],[(DataBa['Mean'].iloc[i] - Quadratic(DataBa['Energy (keV)'].iloc[i],*ParamsLinear)) for i in range(len(DataBa['Energy (keV)']))],fmt='ro',label="Barium-133", markersize=7.5,yerr=DataBa['Error'].values)
plt.errorbar(DataCo['Energy (keV)'],[(DataCo['Mean'].iloc[i] - Quadratic(DataCo['Energy (keV)'].iloc[i],*ParamsLinear)) for i in range(len(DataCo['Energy (keV)']))],fmt='bo',label="Cobalt-60", markersize=7.5,yerr=DataCo['Error'].values)
plt.errorbar(DataNa['Energy (keV)'],[(DataNa['Mean'].iloc[i] - Quadratic(DataNa['Energy (keV)'].iloc[i],*ParamsLinear)) for i in range(len(DataNa['Energy (keV)']))],fmt='yo',label="Sodium-22", markersize=7.5,yerr=DataNa['Error'].values)
plt.title('Energy calibration residual plot',**TitleFont)
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Channel Number Residual',**AxTitleFont)
plt.legend()
plt.show()


ChiStats = ChiSqFunc(list(FitDF['Mean']),list(Quadratic(FitDF['Energy (keV)'],*ParamsLinear)),list(FitDF['Error']),ParamsLinear)
print(ChiStats)





#TODO: Add the energy resolution plot onto this .py file as we have the energies of the peaks as well as the FWHM so we can now produce a nice plot :)



'''







DataCo['Mean'] = [8231,9347]
DataCo['Energy (keV)'] = [1173,1332]

DataNa['Mean'] = [3585,8942]
DataNa['Energy (keV)'] = [511,1274]

DataBa['Mean'] = [570,1132,1574,1948,2135,2509,2705]
DataBa['Energy (keV)'] = [81,161,223,276,303,356,384]


ChannelNos = DataCo['Mean'] + DataNa['Mean'] + DataBa['Mean']
Energies = np.array(DataCo['Energy (keV)'] + DataNa['Energy (keV)'] + DataBa['Energy (keV)'])

Params, Errors = sciopt.curve_fit(linear,Energies,ChannelNos)
print(Params)

rough_EC_a = Params[0]
rough_EC_b = Params[1]

Energies = (df['channel number'] - rough_EC_b)/(rough_EC_a) # keV

plt.semilogy(Energies,df['count number'])
TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'22'}
plt.xlabel('Gamma Energy keV',**AxTitleFont)
plt.ylabel('Log-10 of count number',**AxTitleFont)
plt.title('Mystery Source 2 hour run',**TitleFont)
plt.xlim(0,2270)
plt.show()

plt.plot(DataBa['Energy (keV)'],DataBa['Mean'], 'ro',label="Barium-133", markersize=7.5)
plt.plot(DataCo['Energy (keV)'],DataCo['Mean'],'bo',label="Cobalt-60", markersize=7.5)
plt.plot(DataNa['Energy (keV)'],DataNa['Mean'],'yo',label="Sodium-22", markersize=7.5)

plt.plot(Energies,linear(Energies,*Params))
plt.title('Energy calibration plot with fit',**TitleFont)
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Channel Number',**AxTitleFont)
plt.legend()
plt.show()

Residuals = []

for i in range(len(ChannelNos)):
    Residuals.append(ChannelNos.iloc[i] - linear(Energies.iloc[i],*Params))

plt.plot(DataBa['Energy (keV)'],[(DataBa['Mean'].iloc[i] - linear(DataBa['Energy (keV)'].iloc[i],*Params)) for i in range(len(DataBa['Energy (keV)']))],'ro',label="Barium-133", markersize=7.5)
plt.plot(DataCo['Energy (keV)'],[(DataCo['Mean'].iloc[i] - linear(DataCo['Energy (keV)'].iloc[i],*Params)) for i in range(len(DataCo['Energy (keV)']))],'bo',label="Cobalt-60", markersize=7.5)
plt.plot(DataNa['Energy (keV)'],[(DataNa['Mean'].iloc[i] - linear(DataNa['Energy (keV)'].iloc[i],*Params)) for i in range(len(DataNa['Energy (keV)']))],'yo',label="Sodium-22", markersize=7.5)
plt.title('Energy calibration residual plot',**TitleFont)
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Channel Number Residual',**AxTitleFont)
plt.legend()
plt.show()
'''
