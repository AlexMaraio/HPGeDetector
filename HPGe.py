import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
import scipy as sc
import pandas as pd

#from lmfit import Model, Parameters

'''
 - Produce code to perform energy calibration (fit energy peak values from literature against channel number, linear fit)
 - Fit constrained gaussian profiles to selected data to extract mean channel number and standard deviation
 - Use linear fit parameters to convert channel numbers to energies
 - produce nice plots, correctly labelled.
'''

'''
First thing we need is to accept a csv from the MCA containing values of count vs channel number
To certain parts of this data we should fit a gaussian, extracting the mean value and imposing that the mean channel number corresponds to the energy of the peak.
This should be done for as many peaks, from as many known sources as possible to inform the linear fit as best we can.
We need to find a source from literature for the energy peak values.
'''

def gauss(x,amplitude,mean,sigma,offset):
    return amplitude*np.exp(-((x-mean)**2.0)/(2.0*(sigma**2.0))) + offset
    
BgNoShield = "WeekendNoShield_001_eh_1"
BgNoShieldDF = pd.read_csv(BgNoShield+ ".dat", sep=r"\s+",names = ['channel number','count number'])

MysterySource = "MysterySource-24HrRun_001_eh_1"
df2 = pd.read_csv(MysterySource+ ".dat", sep=r"\s+",names = ['channel number','count number'])


a = 7.05159811
b = -1.02288025

BgNoShieldDF['Energies'] = ( BgNoShieldDF['channel number'] - b )/a 


#print('Total number of counts: ' + str(sum(df['count number'])))

#plt.semilogy(df2['channel number'],df2['count number'])

plt.semilogy(BgNoShieldDF['Energies'],BgNoShieldDF['count number'])
TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'22'}
plt.xlabel('Channel Number',**AxTitleFont)
plt.ylabel('Log-10 of count number',**AxTitleFont)
plt.title('Weekend Background Without Shield',**TitleFont)
plt.xlim(0)
plt.show()








plt.semilogy(BgNoShieldDF['Energies'],BgNoShieldDF['count number'])
TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'22'}
plt.xlabel('Channel Number',**AxTitleFont)
plt.ylabel('Log-10 of count number',**AxTitleFont)
plt.title('Weekend Background Without Shield',**TitleFont)
plt.xlim(0)
plt.show()
