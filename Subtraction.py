import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
import scipy as sc
import pandas as pd

'''
First thing we need is to accept a csv from the MCA containing values of count vs channel number
To certain parts of this data we should fit a gaussian, extracting the mean value and imposing that the mean channel number corresponds to the energy of the peak.
This should be done for as many peaks, from as many known sources as possible to inform the linear fit as best we can.
We need to find a source from literature for the energy peak values.
'''


BgWithShield = "WeekendBgWithShield_003_eh_1"
BgWithShieldDF = pd.read_csv(BgWithShield+ ".dat", sep=r"\s+",names = ['channel number','count number'])

BgWithShieldDFCobalt = BgWithShieldDF.copy()

'''
MysterySource = "MysterySource-24HrRun_001_eh_1"
MysterySourceDF2 = pd.read_csv(MysterySource+ ".dat", sep=r"\s+",names = ['channel number','count number'])
'''

Sodium = 'Sodium22-24HrData_002_eh_1'
SodiumDF = pd.read_csv(Sodium+ ".dat", sep=r"\s+",names = ['channel number','count number'])

Barium = 'Barium133-24HrRun_006_eh_1'
BariumDF = pd.read_csv(Barium+ ".dat", sep=r"\s+",names = ['channel number','count number'])

Cobalt = 'Cobalt60-2HrRun_007_eh_1'
CobaltDF = pd.read_csv(Cobalt+ ".dat", sep=r"\s+",names = ['channel number','count number'])

# BgWithShieldDF['count number'] = BgWithShieldDF['count number'] / ( 239068  ) * 86400

BgWithShieldDF['count number'] = BgWithShieldDF['count number'] / ( 239068  ) * 86400
BgWithShieldDFCobalt['count number'] = BgWithShieldDF['count number'] / ( 239068  ) * (86400/12.0)

#MysterySourceDF2['count number'] = MysterySourceDF2['count number'] - BgWithShieldDF['count number']

SodiumDF['count number'] = SodiumDF['count number'] - BgWithShieldDF['count number']
BariumDF['count number'] = BariumDF['count number'] - BgWithShieldDF['count number']
CobaltDF['count number'] = CobaltDF['count number'] - BgWithShieldDFCobalt['count number']

a = 7.01164642
b = 8.1697501

#MysterySourceDF2['Energies'] = ( MysterySourceDF2['channel number'] - b )/a

SodiumDF['Energies'] = ( SodiumDF['channel number'] - b )/a
BariumDF['Energies'] = ( BariumDF['channel number'] - b )/a
CobaltDF['Energies'] = ( CobaltDF['channel number'] - b )/a

plt.figure(1)
plt.semilogy(SodiumDF['Energies'],SodiumDF['count number'],label="Sodium 22 data")
TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'22'}
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Log-10 of count number',**AxTitleFont)
plt.title('Sodium 22, no background, energies',**TitleFont)
plt.xlim(0)
plt.ylim(1)
plt.legend()


plt.figure(2)
plt.semilogy(BariumDF['Energies'],BariumDF['count number'],label="Barium 133 data")
TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'22'}
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Log-10 of count number',**AxTitleFont)
plt.title('Barium 133, no background, energies',**TitleFont)
plt.xlim(0)
plt.ylim(1)
plt.legend()


plt.figure(3)
plt.semilogy(CobaltDF['Energies'],CobaltDF['count number'],label="Cobalt 60 data")
TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'22'}
plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
plt.ylabel('Log-10 of count number',**AxTitleFont)
plt.title('Cobalt 60, no background, energies',**TitleFont)
plt.xlim(0)
plt.ylim(1)
plt.legend()
plt.show()
#plt.semilogy(BgWithShieldDF['channel number'],BgWithShieldDF['count number'],label="With Shield")
#plt.semilogy(BgWithoutShieldDF2['channel number'],BgWithoutShieldDF2['count number'],label="Without Shield2")
