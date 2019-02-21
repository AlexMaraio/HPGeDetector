import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sciopt
import seaborn as sns
import pandas as pd
sns.set()

source = "MysterySource-2HrRun_001_eh_1"
df = pd.read_csv(source+ ".dat", sep=r"\s+",names = ['channel number','count number'])

TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'22'}


def Linear(E,b,c):
    return E*b + c

def Quadratic(E,a,b,c):
    return E*a**2.0 + E*b + c


FitDF = pd.read_csv('Gamma_Peak_Stats_and_Params.csv', names=['Element', 'Fit type', 'Peak number',	'Peak type', 'Energy (keV)', 'Resolution', 'Min', 'Max', 'Mean', 'A', 'Sigma', 'Error' ,'a', 'b', 'chisq', 'Reduced chisq', 'Probchisq'],usecols=list(range(1,18)))

FitDF = FitDF[FitDF['Fit type'] == 'Linear']
FitDF = FitDF[FitDF['Peak type'] == 'Zoom']

FitDF['Mean'] = FitDF['Mean'].astype(float)

FitDF['Energy (keV)'] = FitDF['Energy (keV)'].astype(float)

plt.plot(FitDF['Energy (keV)'],FitDF['Mean'],'bo')

plt.show()

DataBa = FitDF[FitDF['Element'] == 'Barium133']
DataCo = FitDF[FitDF['Element'] == 'Cobalt60']
DataNa = FitDF[FitDF['Element'] == 'Sodium22']



ParamsLinear, ErrorsLinear = sciopt.curve_fit(Linear,FitDF['Energy (keV)'],FitDF['Mean'])
print(ParamsLinear)

rough_EC_a = ParamsLinear[0]
rough_EC_b = ParamsLinear[1]

Energies = (FitDF['Mean'] - rough_EC_b)/(rough_EC_a) # keV

plt.plot(DataBa['Energy (keV)'],DataBa['Mean'], 'ro',label="Barium-133", markersize=10)
plt.plot(DataCo['Energy (keV)'],DataCo['Mean'],'bo',label="Cobalt-60", markersize=10)
plt.plot(DataNa['Energy (keV)'],DataNa['Mean'],'yo',label="Sodium-22", markersize=10)

plt.plot(Energies,Linear(Energies,*ParamsLinear),label='Fit')
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
plt.plot(DataBa['Energy (keV)'],[(DataBa['Mean'].iloc[i] - Linear(DataBa['Energy (keV)'].iloc[i],*ParamsLinear)) for i in range(len(DataBa['Energy (keV)']))],'ro',label="Barium-133", markersize=10)
plt.plot(DataCo['Energy (keV)'],[(DataCo['Mean'].iloc[i] - Linear(DataCo['Energy (keV)'].iloc[i],*ParamsLinear)) for i in range(len(DataCo['Energy (keV)']))],'bo',label="Cobalt-60", markersize=10)
plt.plot(DataNa['Energy (keV)'],[(DataNa['Mean'].iloc[i] - Linear(DataNa['Energy (keV)'].iloc[i],*ParamsLinear)) for i in range(len(DataNa['Energy (keV)']))],'yo',label="Sodium-22", markersize=10)
plt.title('Energy calibration residual plot',**TitleFont)
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Channel Number Residual',**AxTitleFont)
plt.legend()
plt.show()








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

plt.plot(DataBa['Energy (keV)'],DataBa['Mean'], 'ro',label="Barium-133", markersize=10)
plt.plot(DataCo['Energy (keV)'],DataCo['Mean'],'bo',label="Cobalt-60", markersize=10)
plt.plot(DataNa['Energy (keV)'],DataNa['Mean'],'yo',label="Sodium-22", markersize=10)

plt.plot(Energies,linear(Energies,*Params))
plt.title('Energy calibration plot with fit',**TitleFont)
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Channel Number',**AxTitleFont)
plt.legend()
plt.show()

Residuals = []

for i in range(len(ChannelNos)):
    Residuals.append(ChannelNos.iloc[i] - linear(Energies.iloc[i],*Params))

plt.plot(DataBa['Energy (keV)'],[(DataBa['Mean'].iloc[i] - linear(DataBa['Energy (keV)'].iloc[i],*Params)) for i in range(len(DataBa['Energy (keV)']))],'ro',label="Barium-133", markersize=10)
plt.plot(DataCo['Energy (keV)'],[(DataCo['Mean'].iloc[i] - linear(DataCo['Energy (keV)'].iloc[i],*Params)) for i in range(len(DataCo['Energy (keV)']))],'bo',label="Cobalt-60", markersize=10)
plt.plot(DataNa['Energy (keV)'],[(DataNa['Mean'].iloc[i] - linear(DataNa['Energy (keV)'].iloc[i],*Params)) for i in range(len(DataNa['Energy (keV)']))],'yo',label="Sodium-22", markersize=10)
plt.title('Energy calibration residual plot',**TitleFont)
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Channel Number Residual',**AxTitleFont)
plt.legend()
plt.show()
'''
