import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sciopt
import seaborn as sns
import pandas as pd
sns.set()

source = "MysterySource-2HrRun_001_eh_1"
df = pd.read_table(source+ ".dat", sep=r"\s+",names = ['channel number','count number'])

TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'22'}


def linear(E,a,b):
    return E*a + b


ChannelCo = [8231,9347]
EnergyCo = [1173,1332]

ChannelNa = [3585,8942]
EnergyNa = [511,1274]

ChannelBa = [570,1132,1574,1948,2135,2509,2705]
EnergyBa = [81,161,223,276,303,356,384]


ChannelNos = ChannelCo + ChannelNa + ChannelBa
Energies = np.array(EnergyCo + EnergyNa + EnergyBa)

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

plt.plot(EnergyBa,ChannelBa, 'ro',label="Barium-133", markersize=10)
plt.plot(EnergyCo,ChannelCo,'bo',label="Cobalt-60", markersize=10)
plt.plot(EnergyNa,ChannelNa,'yo',label="Sodium-22", markersize=10)

plt.plot(Energies,linear(Energies,*Params))
plt.title('Energy calibration plot with fit',**TitleFont)
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Channel Number',**AxTitleFont)
plt.legend()
plt.show()

Residuals = []

for i in range(len(ChannelNos)):
    Residuals.append(ChannelNos[i] - linear(Energies[i],*Params))

plt.plot(EnergyBa,[(ChannelBa[i] - linear(EnergyBa[i],*Params)) for i in range(len(EnergyBa))],'ro',label="Barium-133", markersize=10)
plt.plot(EnergyCo,[(ChannelCo[i] - linear(EnergyCo[i],*Params)) for i in range(len(EnergyCo))],'bo',label="Cobalt-60", markersize=10)
plt.plot(EnergyNa,[(ChannelNa[i] - linear(EnergyNa[i],*Params)) for i in range(len(EnergyNa))],'yo',label="Sodium-22", markersize=10)
plt.title('Energy calibration residual plot',**TitleFont)
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Channel Number Residual',**AxTitleFont)
plt.legend()
plt.show()

