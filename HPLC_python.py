import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy
from peakutils import baseline
from scipy.signal import find_peaks
#from scipy.optimize import curve_fit
#from scipy.stats import exponnorm
plt.rcParams.update({'font.size': 16})

############################## DEFINING RATE OF CHANGE ###################################

def rate_change(c1,c0):
	return np.sum(np.abs(c1-c0)/np.abs((np.max([c1,c0],axis=0)+0.000000001)))/((dt*counting))

def rate_change_error(c1,c0,c1_error,c0_error,counting):
	sumeta=0
	for i in range(len(c1)):
		if np.max([c1[i],c0[i]],axis=0)==c1[i]:
			maximet = c1_error[i]
		else:
			maximet = c0_error[i]
		RC0 = (np.abs(c1[i]-c0[i])/(np.max([c1[i],c0[i]],axis=0)+0.000000001))/(dt*counting)
		sumeta += (RC0*(np.sqrt( (c0_error[i]**2+c1_error[i]**2)/((c1[i]-c0[i])**2+0.00000001) + (maximet/(np.max([c1[i],c0[i]],axis=0)+0.0000001)) )))**2
	return np.sqrt(sumeta)

################################## PARAMETERS & VARIABLES ################################
first_half_name = 'CE_rep'
second_half_name = '_C1_to_C15.xlsx'
#first_half_name = 'Moran_CE_long_C3_C45_rep'
#second_half_name = '.xlsx'
num_replicates = 5

dt = 1 # Measurements every 'dt' cycles, if you measure every cycle dt=1
prominence_peaks = 5. # Required hight of peaks with respect to surroundings, there are other parameters to select peaks, but this is the most important one
#Problem: for very little 'prominence_peaks' the identification & tracking seems to fail. Above 3 seems ok.


peaks_removed = [[],[],[],[],[],[16,17]] # Run the code first, identify the number of peaks you want to remove and type it here
#peaks_removed = [[],[],[],[],[],[]]
show_peak_number = True
cycle_plot = 12	# Indicate which cycle do you want to plot the absorbance vs RT with mean and error
cycle_marked_black = None # Indicate if you want to highlight a cycle among all, put 'None' otherwise
replicate_plot = num_replicates # Index of the replicate you want to plot. If you wan to plot the mean, put 'replicate_plot = num_replicates'

baseline_cut_location = 2 #location of the cut, 2 means we take away the first half, 3 means take away the first third...
baseline_polinomial_order = 3 #order of the polinomial for the baseline correction
max_distance = 40 # Maximal distance for peaks to be considered the same. If within cycles peaks are displaced in the RT axis more than this number they will be considered different. Two peaks that are different and are separated less than this distance may be considered the same if one disappears.

RT_resolution = 0.006667

################################## READING DATA ##########################################

#### Read data
df = pd.read_excel(first_half_name+str(1)+second_half_name)
num_RT = len(df)
num_cycles = len(df.columns)-2

C,RT = [],[]
for k in range(num_RT):
	RT.append(df.values[k][1])
for i in range(num_replicates):
	df = pd.read_excel(first_half_name+str(i+1)+second_half_name)
	C.append([])
	for j in range(num_cycles):
		C[i].append([])
		for k in range(num_RT):
			C[i][j].append(df.values[k][j+2])
df = 0

print ('# Replicates = ',len(C))
print ('# Cycles = ',len(C[0]))
print ('# RT = ',len(C[0][0]))

#### Calculate mean
C.append([])
i = len(C)-1
for j in range(num_cycles):
	C[i].append([])
	for k in range(num_RT):
		sumeta = 0
		for ii in range(num_replicates):
			sumeta += C[ii][j][k]/num_replicates
		C[i][j].append(sumeta)
		
#### Calculate error
C.append([])
i = len(C)-1
for j in range(num_cycles):
	C[i].append([])
	for k in range(num_RT):
		sumeta = 0
		for ii in range(num_replicates):
			sumeta += ((C[ii][j][k]-C[i-1][j][k])**2)/(num_replicates-1)
		C[i][j].append(2.*np.sqrt(sumeta))
					
################################# DATA MODIFICATIONS #####################################

### Save original data
C0 = []
for i in range(len(C)):
	C0.append([])
	for j in range(len(C[i])):
		C0[i].append(C[i][j])

### Data cut
for i in range(len(C)):
	for j in range(len(C[i])):
		C[i][j] = C[i][j][int(len(C[i][j])/baseline_cut_location):]

# Base line correction
C = np.array(C)
for i in range(len(C)):
	for j in range(len(C[i])):
		bl = baseline(C[i][j],deg=baseline_polinomial_order)
		C[i][j] += -bl
C = C.tolist()

### Data cut back
for i in range(len(C)):
	for j in range(len(C[i])):
		C[i][j] = np.concatenate((C0[i][j][:int(len(C0[i][j])/baseline_cut_location)],C[i][j]),axis=0)

### Calculate new mean
for i in range(num_cycles):
	sumeta_1=0
	for j in range(num_replicates):
		sumeta_1 += C[j][i]/num_replicates
	C[num_replicates][i] = sumeta_1
	
### Calculate new error
for i in range(num_cycles):
	sumeta_2=0
	for j in range(num_replicates):
		sumeta_2 += (1/(num_replicates-1))*(C[j][i]-C[num_replicates][i])**2
	C[num_replicates+1][i] = 2.*np.sqrt(sumeta_2)

################################## PLOT DATA #############################################
### Plot original data, for a particular cycle, each replicate, the mean and error at 95%
fig, ax = plt.subplots()
C0 = np.array(C0)
for i in range(num_replicates):
	if i==0:
		plt.plot(RT,C0[i][cycle_plot],color='b',label='Original')
	else:
		plt.plot(RT,C0[i][cycle_plot],color='b')#,alpha=1/float(i+1))
plt.plot(RT,C0[num_replicates][cycle_plot],color='black')
#ax.errorbar(RT,C0[num_replicates][cycle_plot], yerr=C0[num_replicates+1][cycle_plot], fmt='-o',c='black')
plt.fill_between(RT,C0[num_replicates][cycle_plot]-C0[num_replicates+1][cycle_plot],C0[num_replicates][cycle_plot]+C0[num_replicates+1][cycle_plot],facecolor='b',alpha=0.2)
### Plot modified data, for a particular cycle, each replicate, the mean and error at 95%
C = np.array(C)		
for i in range(num_replicates):
	if i==0:
		plt.plot(RT,C[i][cycle_plot]+50,color='r',label='Correction')
	else:
		plt.plot(RT,C[i][cycle_plot]+50,color='r')#,alpha=1/float(i+1))
plt.plot(RT,C[num_replicates][cycle_plot]+50,color='black')
plt.fill_between(RT,50+C[num_replicates][cycle_plot]-C[num_replicates+1][cycle_plot],50+C[num_replicates][cycle_plot]+C[num_replicates+1][cycle_plot],facecolor='r',alpha=0.2)
plt.legend()
plt.title('# Cycle = '+str(cycle_plot))
plt.show()

### Plot all data as HeatMap
C_map = []
eps = -C[num_replicates].min()+1.
for i in range(num_cycles):
	C_map.append([])
	for j in range(num_RT):
		C_map[i].append(np.log(C[num_replicates][i][j]+eps))
plt.imshow(C_map, cmap='Reds', aspect='auto')
plt.xlabel('RT')
plt.ylabel('Cycle')
plt.title('log (Absorbance)')
plt.colorbar()
plt.show()

################################# PICK FINDER & INTEGRATION ##############################
#height: required hight of picks, prominence: required hight with respect to others, width: required width, threshold: required distance from neighbour pick

num_color = 0
peaks_index,peaks_integral,peaks_integral_error = [],[],[]

for i in range(num_replicates + 1):
	peaks_index.append([])
	peaks_integral.append([])
	for j in range(num_cycles):
		peaks, info = find_peaks(C[i][j],height=0,prominence=prominence_peaks,width=0.,threshold=0.)
		peaks_index[i].append(peaks)
		peaks_integral[i].append([])
		#lower = info['left_bases'] #some peaks exceed their limit, so is better let_ips
		#upper = info['right_bases']
		lower = info['left_ips']
		upper = info['right_ips']
		if i == num_replicates:
			peaks_integral_error.append([])

		for k in range(len(peaks_index[i][j])):
			#peaks_integral[i][j].append(0.006667*np.trapz(C[i][j][int(lower[k]):int(upper[k])])) # automatic tepezoid integration
			manual_integral,manual_integral_error = 0,0
			#for m in range(peaks_index[i][j][k]-1,peaks_index[i][j][k]+0): # loop without integration
			for m in range(int(lower[k]),int(upper[k])):
				manual_integral += (RT_resolution*(C[i][j][m]+C[i][j][m+1])/2)
				if i == num_replicates:
					manual_integral_error += ((RT_resolution*(C[num_replicates+1][j][m]+C[num_replicates+1][j][m+1])/2))**2
			peaks_integral[i][j].append(manual_integral)
			if i == num_replicates:
				manual_integral_error = np.sqrt(manual_integral_error)
				peaks_integral_error[j].append(manual_integral_error)
		
################################## PLOT PEAKS INTEGRAL ###################################
			if i==replicate_plot and j==cycle_plot:
				if num_color == 0:
					plt.plot(C[i][j],c='black',zorder=0)
					plt.fill_between(range(int(lower[k]),int(upper[k])),0,C[i][j][int(lower[k]):int(upper[k])],color='C'+str(num_color))
				plt.scatter(np.array(RT[peaks_index[i][j][k]])/RT_resolution,C[i][j][peaks_index[i][j][k]],c='C'+str(num_color),zorder=1)
				plt.fill_between(range(int(lower[k]),int(upper[k])),0,C[i][j][int(lower[k]):int(upper[k])],color='C'+str(num_color))
				num_color += 1
				plt.title('# Cycle = '+str(cycle_plot))
		plt.show()

########### PEACK IDENTIFICATOR & ADDITION OF ZEROS FOR PROPER COMPARISON ################

### Add zero values to later cycles for those peaks appearing in initial cycles
for sample in range(num_replicates+1):
	for cycle in range(num_cycles-1):
		
		#Adding zeros to peaks missing from cycle 'i' to 'i+1'
		for i in peaks_index[sample][cycle]:
			minimum_1 = 1000
			for j in peaks_index[sample][cycle+1]:
				minimum_2 = 1000
				if np.abs(i-j)<minimum_1:
					minimum_1 = np.abs(i-j)
					jj=j
			for k in peaks_index[sample][cycle]:
				if np.abs(jj-k)<minimum_2 and np.abs(jj-k)<max_distance:
					minimum_2 = np.abs(jj-k)
					kk=k
			if kk==i:
				pass
			else:
				peaks_index[sample][cycle+1]= np.append(peaks_index[sample][cycle+1],[i],axis=0)
				peaks_integral[sample][cycle+1]= np.append(peaks_integral[sample][cycle+1],[0],axis=0)
				peaks_integral_error[cycle+1]= np.append(peaks_integral_error[cycle+1],[0],axis=0)	
		
		#Ordering from min RT to max RT
		for n in range(len(peaks_index[sample][cycle+1])):
			for m in range(0, len(peaks_index[sample][cycle+1])-n-1):
				if peaks_index[sample][cycle+1][m] > peaks_index[sample][cycle+1][m+1]:
					peaks_index[sample][cycle+1][m], peaks_index[sample][cycle+1][m+1] = peaks_index[sample][cycle+1][m+1], peaks_index[sample][cycle+1][m]
					peaks_integral[sample][cycle+1][m], peaks_integral[sample][cycle+1][m+1] = peaks_integral[sample][cycle+1][m+1], peaks_integral[sample][cycle+1][m]
					if sample == num_replicates:
						peaks_integral_error[cycle+1][m], peaks_integral_error[cycle+1][m+1] = peaks_integral_error[cycle+1][m+1], peaks_integral_error[cycle+1][m]
		
		#Adding zeros to peaks missing from cycle 'i' to 'i+1'
		for i in peaks_index[sample][cycle+1]:
			minimum_1 = 1000
			for j in peaks_index[sample][cycle]:
				minimum_2 = 1000
				if np.abs(i-j)<minimum_1:
					minimum_1 = np.abs(i-j)
					jj=j	
			for k in peaks_index[sample][cycle+1]:
				if np.abs(jj-k)<minimum_2 and np.abs(jj-k)<max_distance:
					minimum_2 = np.abs(jj-k)
					kk=k
			if kk==i:
				pass
			else:
				peaks_index[sample][cycle]= np.append(peaks_index[sample][cycle],[i],axis=0)
				peaks_integral[sample][cycle]= np.append(peaks_integral[sample][cycle],[0],axis=0)
				if sample == num_replicates:
					peaks_integral_error[cycle]= np.append(peaks_integral_error[cycle],[0],axis=0)
		
		#Ordering again
		for n in range(len(peaks_index[sample][cycle])):
			for m in range(0, len(peaks_index[sample][cycle])-n-1):
				if peaks_index[sample][cycle][m] > peaks_index[sample][cycle][m+1]:
					peaks_index[sample][cycle][m], peaks_index[sample][cycle][m+1] = peaks_index[sample][cycle][m+1], peaks_index[sample][cycle][m]
					peaks_integral[sample][cycle][m], peaks_integral[sample][cycle][m+1] = peaks_integral[sample][cycle][m+1], peaks_integral[sample][cycle][m]	
					if sample == num_replicates:
						peaks_integral_error[cycle][m], peaks_integral_error[cycle][m+1] = peaks_integral_error[cycle][m+1], peaks_integral_error[cycle][m]	

### Repetition backwards to add zero values to initial cycles for those peaks appearing in the last cycles
for sample in range(num_replicates+1):
	for cycle in range(num_cycles-2,-1,-1):

		for i in peaks_index[sample][cycle]:
			minimum_1 = 1000
			for j in peaks_index[sample][cycle+1]:
				minimum_2 = 1000
				if np.abs(i-j)<minimum_1:
					minimum_1 = np.abs(i-j)
					jj=j
			for k in peaks_index[sample][cycle]:
				if np.abs(jj-k)<minimum_2 and np.abs(jj-k)<max_distance:
					minimum_2 = np.abs(jj-k)
					kk=k
			if kk==i:
				pass
			else:
				peaks_index[sample][cycle+1]= np.append(peaks_index[sample][cycle+1],[i],axis=0)
				peaks_integral[sample][cycle+1]= np.append(peaks_integral[sample][cycle+1],[0],axis=0)
				if sample == num_replicates:
					peaks_integral_error[cycle+1]= np.append(peaks_integral_error[cycle+1],[0],axis=0)	
		
		for n in range(len(peaks_index[sample][cycle+1])):
			for m in range(0, len(peaks_index[sample][cycle+1])-n-1):
				if peaks_index[sample][cycle+1][m] > peaks_index[sample][cycle+1][m+1]:
					peaks_index[sample][cycle+1][m], peaks_index[sample][cycle+1][m+1] = peaks_index[sample][cycle+1][m+1], peaks_index[sample][cycle+1][m]
					peaks_integral[sample][cycle+1][m], peaks_integral[sample][cycle+1][m+1] = peaks_integral[sample][cycle+1][m+1], peaks_integral[sample][cycle+1][m]
					if sample == num_replicates:
						peaks_integral_error[cycle+1][m], peaks_integral_error[cycle+1][m+1] = peaks_integral_error[cycle+1][m+1], peaks_integral_error[cycle+1][m]
		
		for i in peaks_index[sample][cycle+1]:
			minimum_1 = 1000
			for j in peaks_index[sample][cycle]:
				minimum_2 = 1000
				if np.abs(i-j)<minimum_1:
					minimum_1 = np.abs(i-j)
					jj=j	
			for k in peaks_index[sample][cycle+1]:
				if np.abs(jj-k)<minimum_2 and np.abs(jj-k)<max_distance:
					minimum_2 = np.abs(jj-k)
					kk=k
			if kk==i:
				pass
			else:
				peaks_index[sample][cycle]= np.append(peaks_index[sample][cycle],[i],axis=0)
				peaks_integral[sample][cycle]= np.append(peaks_integral[sample][cycle],[0],axis=0)
				if sample == num_replicates:
					peaks_integral_error[cycle]= np.append(peaks_integral_error[cycle],[0],axis=0)
		
		for n in range(len(peaks_index[sample][cycle])):
			for m in range(0, len(peaks_index[sample][cycle])-n-1):
				if peaks_index[sample][cycle][m] > peaks_index[sample][cycle][m+1]:
					peaks_index[sample][cycle][m], peaks_index[sample][cycle][m+1] = peaks_index[sample][cycle][m+1], peaks_index[sample][cycle][m]
					peaks_integral[sample][cycle][m], peaks_integral[sample][cycle][m+1] = peaks_integral[sample][cycle][m+1], peaks_integral[sample][cycle][m]	
					if sample == num_replicates:
						peaks_integral_error[cycle][m], peaks_integral_error[cycle][m+1] = peaks_integral_error[cycle][m+1], peaks_integral_error[cycle][m]	

############################# PLOT PEAKS IDENTIFIED ######################################

fig, ax = plt.subplots()
colors = plt.cm.winter(np.linspace(0,1,num_cycles))
peak_number = 0
#colors = plt.cm.hsv(np.linspace(0,1,16))
for i in range(num_cycles):
	if i==cycle_marked_black:
		plt.plot(RT,C[replicate_plot][i]+0*i,c='black',zorder=0)
	else:	
		plt.plot(RT,C[replicate_plot][i]+0*i,c=colors[i],zorder=0)
	num_color = 0
	for j in peaks_index[replicate_plot][i]:
		plt.scatter(RT[j],C[replicate_plot][i][j]+0*i,c='C'+str(num_color),zorder=1)
		if i==num_cycles-1 and show_peak_number==True:
			ax.annotate(str(peak_number), (RT[j],C[replicate_plot][i][j]),fontsize=7)
			peak_number += 1
		num_color += 1
print ('# Peaks identified in replicate plotted = ',len(peaks_index[replicate_plot][i]))
plt.show()

############################### PLOT EVOLUTION ###########################################

### Remove peaks manually identified as wrong, to correct the synchrony plot and the rate of change
peaks_integral_0,peaks_integral_error_0 = [],[]
for i in range(num_replicates+1):
	peaks_integral_0.append([])
	for j in range(num_cycles):
		peaks_integral_0[i].append([])
		if i==num_replicates:
			peaks_integral_error_0.append([])
		for k in range(len(peaks_integral[i][j])):
			if k in peaks_removed[i]:
				pass
			else:
				peaks_integral_0[i][j].append(peaks_integral[i][j][k])
				if i==num_replicates:
					peaks_integral_error_0[j].append(peaks_integral_error[j][k])

peaks_integral,peaks_integral_error = peaks_integral_0,peaks_integral_error_0
		
### Correcting the data for peaks that disappear and appear suddenly due to errors and putting the evolution of each peak in a different vector
peak_evolution,peak_evolution_error,cycle_evolution = [],[],[]
for k in range(num_replicates+1):
	peak_evolution.append([])
	for i in range(len(peaks_integral[k][0])):
		peak_evolution[k].append([])
		if k==num_replicates:
			peak_evolution_error.append([])
		for j in range(num_cycles):
			if j!=0 and j!=(num_cycles-1) and peaks_integral[k][j][i]==0 and peaks_integral[k][j-1][i]!=0 and peaks_integral[k][j+1][i]!=0:	
				peak_evolution[k][i].append((peaks_integral[k][j-1][i]+peaks_integral[k][j+1][i])/2.)
				if k==num_replicates:
					peak_evolution_error[i].append((peaks_integral_error[j-1][i]+peaks_integral_error[j+1][i])/2.)
			else:
				peak_evolution[k][i].append(peaks_integral[k][j][i])
				if k==num_replicates:	
					peak_evolution_error[i].append(peaks_integral_error[j][i])
				
### Synchrony plot, with errorbars for the mean. The error is propagated from the deviation in the absorbance, it is not the error of the deviation between different synchrony plots. Different synchrony plots may differ because there is some thresholds for detection...

cycle_evolution = np.linspace(0,dt*(num_cycles-1),num=num_cycles).tolist()
num_color = 0
num_plots = max(int(len(peak_evolution[replicate_plot])/10),1)
fig, axs = plt.subplots(num_plots+1, sharex=False, sharey=False)
fig.suptitle('Synchrony')
for i in range(min(len(peak_evolution[replicate_plot]),10)):
	for j in range(num_plots):
		num_color = i+j*10
		axs[j].plot(cycle_evolution,peak_evolution[replicate_plot][i+j*10]/np.max(peak_evolution[replicate_plot][i+j*10]),c='C'+str(num_color))
		if replicate_plot==num_replicates:
			axs[j].fill_between(cycle_evolution,peak_evolution[replicate_plot][i+j*10]/np.max(peak_evolution[replicate_plot][i+j*10])-peak_evolution_error[i+j*10]/np.max(peak_evolution[replicate_plot][i+j*10]),peak_evolution[replicate_plot][i+j*10]/np.max(peak_evolution[replicate_plot][i+j*10])+peak_evolution_error[i+j*10]/np.max(peak_evolution[replicate_plot][i+j*10]),facecolor='C'+str(num_color),alpha=0.1)
for i in range(len(peak_evolution[replicate_plot])-(num_plots)*10):
	num_color = i+num_plots*10
	axs[num_plots-1].plot(cycle_evolution,peak_evolution[replicate_plot][num_plots*10+i]/np.max(peak_evolution[replicate_plot][num_plots*10+i]),c='C'+str(num_color))
	if replicate_plot==num_replicates:
		axs[num_plots-1].fill_between(cycle_evolution,peak_evolution[replicate_plot][num_plots*10+i]/np.max(peak_evolution[replicate_plot][num_plots*10+i])-peak_evolution_error[num_plots*10+i]/np.max(peak_evolution[replicate_plot][num_plots*10+i]),peak_evolution[replicate_plot][num_plots*10+i]/np.max(peak_evolution[replicate_plot][num_plots*10+i])+peak_evolution_error[num_plots*10+i]/np.max(peak_evolution[replicate_plot][num_plots*10+i]),facecolor='C'+str(num_color),alpha=0.1)

### Plot all peaks identified after removal of the wrong ones
peak_number = 0
for i in range(num_cycles):
	axs[num_plots].plot(RT,C[replicate_plot][i]+0*i,c=colors[i],zorder=0)
	num_color = 0
	jj=0
	for j in peaks_index[replicate_plot][i]:
		if jj in peaks_removed[replicate_plot]:
			pass
		else:
			axs[num_plots].scatter(RT[j],C[replicate_plot][i][j]+0*i,c='C'+str(num_color),zorder=1)
			if i==num_cycles-1 and show_peak_number==True:
				axs[num_plots].annotate(str(peak_number), (RT[j],C[replicate_plot][i][j]),fontsize=7)
				peak_number += 1
			num_color += 1
		jj+=1
plt.show()

################################# RATE OF CHANGE #########################################

# STILL NEED TO GIVE THE OPTION TO REMOVE SOME OF THE PEAKS
peaks_removed = []

RateChange, RateChangeError = [],[]
for i in range(num_replicates+1):
	RateChange.append([])
	for j in range(num_cycles-1):
		counting = 0
		for k in range(len(peaks_integral[i][j])):
			if peaks_integral[i][j][k]==0 and peaks_integral[i][j+1][k]==0:
				pass
			else:
				counting += 1
		RateChange[i].append(rate_change(np.array(peaks_integral[i][j]),np.array(peaks_integral[i][j+1])))	
		if i==num_replicates:
			RateChangeError.append(rate_change_error(np.array(peaks_integral[i][j]),np.array(peaks_integral[i][j+1]),np.array(peaks_integral_error[j]),np.array(peaks_integral_error[j+1]),counting))

cycle_RateChange = np.linspace(1,dt*(num_cycles-1),num=num_cycles-1).tolist()
fig, ax = plt.subplots()
plt.xlim(0,dt*num_cycles)
plt.ylim(-0.,.9)
plt.axhline(y=0.1, color='black', linestyle='--',linewidth=1.)
ax.errorbar(cycle_RateChange,RateChange[num_replicates], yerr=RateChangeError, fmt='-o',c='black')
plt.show()

fig, ax = plt.subplots()
plt.xlim(0,dt*num_cycles)
plt.ylim(-0.,.9)
plt.axhline(y=0.1, color='black', linestyle='--',linewidth=1.)
for i in range(num_replicates):
	plt.plot(cycle_RateChange,RateChange[i])
ax.errorbar(cycle_RateChange,RateChange[num_replicates], yerr=RateChangeError, fmt='-o',c='black')
plt.show()