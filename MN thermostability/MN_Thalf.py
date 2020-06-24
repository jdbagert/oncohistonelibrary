########################################################################################
##
##	MN_Thalf.py
##	John Bagert - 2020
##
##	This script calculates the Thalf values from mononucleosome (MN) melt experiments,
##	based on a sigmoidal fit of the dimer-melt section of the curve.
##	It uses fluorescent measurement outputs from a ViiA 7 qPCR instrument. 
##
##	NOTE - The data must be sorted externally by cycle (this script does not sort)
##
########################################################################################

### imports modules, defines functions, other housekeeping
####################################################################

# Imports modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import seaborn as sns
from matplotlib import rc, rcParams
from math import ceil
import pprint
from scipy.optimize import curve_fit

# Housekeeping: Resets variables, closes plots, adjusts figure cosmetics
if 1:
	from IPython import get_ipython
	get_ipython().magic('reset -sf') 
	### Housekeeping
	plt.close('all')
	font = {'family' : 'sans',
			'weight' : 'normal',
			'size' : 10}
	rc('pdf', fonttype=42)
	rc('font',**font)
	rcParams.update({'font.size': 12,
		'font.sans-serif':'Arial'})
	rcParams['mathtext.fontset'] = 'custom'
	rcParams['mathtext.rm'] = 'Arial'
	rcParams['text.color'] = 'black'
	rcParams['axes.labelcolor'] = 'black'
	sns.set_style('ticks')

# Defines sigmoid function for fitting thermal melts
def sigmoid(x, x0, y0,c, k):
     y = y0+c*(1 / (1 + np.exp(-k*(x-x0))))
     return y


### script parameters
####################################################################
infile = 'data.txt' # NOTE! DATA MUST BE SORTED
insamples = 'samples.txt' # defines what wells correspond to what samples
Tvalues = range(25,96) # melt curve from 25-95 C
showplots = 1	# whether or not to show plots
savefig = 1		# whether or not to save figures
replicates = {}	# stores replicate information, builds from insamples file automatically


### reading / storing / normalizing data
#####################################################################

# Initiates variables to hold sample names, replicates info, and sample data
well2samp = {}
fsamps = open(insamples)
fsampsh = fsamps.readline().lower().strip().split('\t')
for l in fsamps:
	d = dict(zip(fsampsh,l.strip().split('\t')))
	well2samp[d['well']] = d['sample']
	if not d['sample'] in replicates:
		replicates[d['sample']] = []
	replicates[d['sample']].append(d['well'])
data = {}
for l in well2samp:	data[l]=[]

# Reads data file and stores well / fluorescence information
f = open(infile)
fh = f.readline().lower().strip().split('\t')
for l in f:
	d = dict(zip(fh,l.strip().split('\t')))
	_wellpos = d['well position']
	data[_wellpos.strip()].append(float(d['tamra'].replace('"','').replace(',','')))

# Normalizes data per well
data_norm = {}
for l in data:
	_max = max(data[l][-30:])
	_mean = np.mean(data[l][30:36])
	data_norm[l] = [(x-_mean)/(_max-_mean) for x in data[l]]	# Normalizes max value and average baseline min value


### Performs sigmoid fits and plots data
####################################################################
if 1:
	for l in replicates:
		_title = l
		_samps = replicates[l]
		plt.figure(facecolor='w',figsize=(4,3))
		Tm = []
		for m in _samps:
			datatofit_y = data_norm[m][20:-13]
			datatofit_x = Tvalues[20:20+len(datatofit_y)]
			ydata = np.array(datatofit_y)
			xdata = np.array(datatofit_x)
			plt.plot(Tvalues,data_norm[m], 'o', ms = 5, mew=.5, mfc = 'none') # plots data to fit
			
			### Fits data
			try: popt, pcov = curve_fit(sigmoid, xdata, ydata,p0=[73,1,1,0])
			except RuntimeError:
				print 'Could not fit %s' %(l)
				continue
			x = np.linspace(xdata[0],xdata[-1],50)
			y = sigmoid(x, *popt)
			
			### Plots fits
			plt.plot(x,y)
			plt.plot([popt[0]]*2,[0,1],'--k',alpha=.33)
			plt.plot([popt[3]]*2,[0,1],'--k',alpha=.33)
			Tm.append(popt[0])
			
		plt.text(43,1,l)
		plt.ylabel('SYPRO Orange (RFU)')
		plt.xlabel('Temperature ($\degree$C)')
		plt.xlim([40,98])
		plt.xticks([40,50,60,70,80,90])
		plt.ylim([-.1,1.1])		
		_Tmavg = np.mean(Tm)
		_Tmstd = np.std(Tm)
		plt.text(_Tmavg,1.01,'%.2f $\pm$ %.2f' %(_Tmavg,_Tmstd),horizontalalignment='center')
		print _title + ': Tm = %.2f +/- %.2f' %(_Tmavg,_Tmstd)
		
		# Puts tickmarks on outside
		ax = plt.gca()
		ax.get_yaxis().set_tick_params(direction='out')
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		# Only show ticks on the left and bottom spines
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')
		# Change axis and tick linewidth
		ax.spines['left'].set_linewidth(1)
		ax.spines['bottom'].set_linewidth(1)
		#ax.spines['bottom'].set_linewidth(2)
		ax.yaxis.set_tick_params(width=1)
		ax.xaxis.set_tick_params(width=1)
		#for label in ax.get_yticklabels():
		#	label.set_horizontalalignment('left')
		ax.tick_params(pad=2,length=4,labelsize=10)
		plt.tight_layout()
	
		if savefig:
			plt.savefig(_title + '.pdf')
		
if showplots: plt.show()
