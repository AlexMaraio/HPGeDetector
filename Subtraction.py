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

# If statements that control what the code does! 
PlotMysterySource = True
PlotKnownSources = True

# Plot controlls
TitleFont = {'size':'23', 'color':'black', 'weight':'bold'} 
AxTitleFont = {'size':'20'}

if PlotMysterySource:
	BgWithShield = "WeekendBgWithShield_003_eh_1"
	BgWithShieldDF = pd.read_csv(BgWithShield+ ".dat", sep=r"\s+",names = ['channel number','count number'])
	
	MysterySource = "MysterySource-24HrRun_001_eh_1"
	MysterySourceDF = pd.read_csv(MysterySource+ ".dat", sep=r"\s+",names = ['channel number','count number'])
	
	# Normalises the timings of the two runs to 24 hours each
	BgWithShieldDF['count number'] = BgWithShieldDF['count number'] / ( 239068  ) * 86400
	
	MysterySourceDF['count number'] = MysterySourceDF['count number'] - BgWithShieldDF['count number']
	
	a = 7.01164642
	b = 8.1697501

	MysterySourceDF['Energies'] = ( BgWithShieldDF['channel number'] - b )/a 
	
	plt.semilogy(MysterySourceDF['Energies'],MysterySourceDF['count number'],label="Mystery source")
	plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
	plt.ylabel('Log-10 of count number',**AxTitleFont)
	plt.title('Mystery Source With Background Subtracted',**TitleFont)
	plt.xlim(0)
	plt.ylim(1)
	plt.legend()
	plt.show()


if PlotKnownSources:
	BgWithShield = "WeekendBgWithShield_003_eh_1"
	BgWithShieldDF = pd.read_csv(BgWithShield+ ".dat", sep=r"\s+",names = ['channel number','count number'])
	
	Sodium = 'Sodium22-24HrData_002_eh_1'
	SodiumDF = pd.read_csv(Sodium+ ".dat", sep=r"\s+",names = ['channel number','count number'])
	
	Barium = 'Barium133-24HrRun_006_eh_1'
	BariumDF = pd.read_csv(Barium+ ".dat", sep=r"\s+",names = ['channel number','count number'])
	
	Cobalt = 'Cobalt60-2HrRun_007_eh_1'
	CobaltDF = pd.read_csv(Cobalt+ ".dat", sep=r"\s+",names = ['channel number','count number'])
	
	# Normalises the background to the correct timings
	BgWithShieldDFCobalt = BgWithShieldDF.copy()
	BgWithShieldDF['count number'] = BgWithShieldDF['count number'] / ( 239068  ) * 86400
	BgWithShieldDFCobalt['count number'] = BgWithShieldDF['count number'] / ( 239068  ) * (86400/12.0)
	
	# Subtracts the background from the original data
	SodiumDF['count number'] = SodiumDF['count number'] - BgWithShieldDF['count number']
	BariumDF['count number'] = BariumDF['count number'] - BgWithShieldDF['count number']
	CobaltDF['count number'] = CobaltDF['count number'] - BgWithShieldDFCobalt['count number']
		
	# Fitting parameters
	a = 7.01164642
	b = 8.1697501	
	
	# Converting channel numbers to energies
	SodiumDF['Energies'] = ( SodiumDF['channel number'] - b )/a
	BariumDF['Energies'] = ( BariumDF['channel number'] - b )/a
	CobaltDF['Energies'] = ( CobaltDF['channel number'] - b )/a
	
	
	# Now plotting the different elements on different log-y plots.
	
	plt.figure(1)
	plt.semilogy(SodiumDF['Energies'],SodiumDF['count number'],label="Sodium 22 data",color='yellow')
	TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
	AxTitleFont = {'size':'22'}
	plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
	plt.ylabel('Log-10 of count number',**AxTitleFont)
	plt.title('Sodium 22, no background, energies',**TitleFont)
	plt.xlim(0)
	plt.ylim(1)
	plt.legend()


	plt.figure(2)
	plt.semilogy(BariumDF['Energies'],BariumDF['count number'],label="Barium 133 data",color='green')
	TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
	AxTitleFont = {'size':'22'}
	plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
	plt.ylabel('Log-10 of count number',**AxTitleFont)
	plt.title('Barium 133, no background, energies',**TitleFont)
	plt.xlim(0)
	plt.ylim(1)
	plt.legend()


	plt.figure(3)
	plt.semilogy(CobaltDF['Energies'],CobaltDF['count number'],label="Cobalt 60 data",color='blue')
	TitleFont = {'size':'25', 'color':'black', 'weight':'bold'} 
	AxTitleFont = {'size':'22'}
	plt.xlabel('Gamma Energy [keV]',**AxTitleFont)
	plt.ylabel('Log-10 of count number',**AxTitleFont)
	plt.title('Cobalt 60, no background, energies',**TitleFont)
	plt.xlim(0)
	plt.ylim(1)
	plt.legend()
	plt.show()

