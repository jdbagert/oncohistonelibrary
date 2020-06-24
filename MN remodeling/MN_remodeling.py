########################################################################################
##
##	MN_remodeling.py
##	John Bagert - 2020
##
##	This script calculates the observed rate values (Kobs) in mononucleosome (MN) 
##	remodeling experiments, based on an exponential decay fit. Kobs values are 
##	written to an output file.
##
##	NOTE - The data must be sorted externally by cycle (this script does not sort)
##
########################################################################################


### imports modules, defines functions, other housekeeping
####################################################################
from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf

### Housekeeping
plt.close('all')
font = {'family' : 'sans',
		'weight' : 'normal',
		'size' : 12}
rc('pdf', fonttype=42)
rc('font',**font)

### Fit functions - exponential decay with a y-int of 1 (normalized data)
def exp_decay(t, kobs):
    f = 1*np.exp(-kobs*t)+yinf
    return f
times = [0,2,4,8,16,32,64] # time-points in remodeling experiments (minutes)
yinf = 0.0 # x-int set to 0


### script parameters, variables
####################################################################

### Input files
fn = 'data.txt' # input data file
f = open(fn)
fh = f.readline().lower().strip().split('\t')

### Whether to save & show figures
savefig = 1 	# whether to save figure
plotfigs = 1	# whether to plot figures
showfig = 1		# whether to show figures

### Defines variables
alldata = {}
kdic = {}
errdic = {}
r2 = {}


### reads & stores file data
####################################################################

### Read file data
for l in f:
	line = l.strip().split('\t')
	d = dict(zip(fh,line))
	_mutant = d['mutant histone']
	_y = line[-1*3*len(times):]
	_y = [float(x) for x in _y]
	alldata[_mutant] = _y
	kdic[_mutant] = []
	r2[_mutant] = []
	
### performs fits 
kguess = [.01]
t = times
t3 = times+times+times
for _mn in alldata:
	y = alldata[_mn]
	popt, pcov = curve_fit(exp_decay,t3,y,kguess)
	perr = np.sqrt(np.diag(pcov))
	kdic[_mn] = popt[0]
	errdic[_mn] = perr
	res = y - exp_decay(t3,popt)
	ss_res = np.sum(res**2)
	ss_tot = np.sum((y-np.mean(y))**2)
	r2[_mn] = 1-(ss_res/ss_tot)


### plots data and outputs Kobs values
####################################################################

### Plots data
if plotfigs:
	for l in alldata:
		_t = np.array(range(0,t[-1],1))
		_fit = exp_decay(_t, kdic[l])

		### Plots figure
		plt.figure(facecolor='w')
		plt.plot(t3,alldata[l],'ok',ms=6,mew=1,mfc='None',label='data')
		plt.xlabel('min')
		plt.ylabel('DNA')
		plt.ylim([0,1.1])
		plt.plot(_t,_fit,'m-',lw=2,label = 'fit')
		plt.legend(loc='upper right')
		plt.title('%s R$^2$ = %.2f' %(l,r2[l]))
		plt.tight_layout()
		if savefig: plt.savefig(l+'.jpg')
		if not showfig: plt.close()

### Writes Kobs values to output file
fout = open('output.txt','w')
fout.write('MN\tkobs\terr\n')
for l in kdic:
	fout.write(l+'\t'+str(kdic[l])+'\t'+str(errdic[l][0])+'\n')
fout.close()

### Show figures
if showfig: plt.show()