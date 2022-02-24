'''
#################################
run_stats.py
################################# 

This model was written by Abygail R. Waggoner(1*) and L. Ilsedore Cleeves(1,2).

This script calculates a flare energy distribution from any number of 
light curves produced by xray_generator.py

(1) University of Virginia, Department of Chemistry
(2) University of virginia, Department of Astronomy
* contact information: arw6qz@virginia.edu

##################################
To use this script: 

1.) Generate light curves with xray_generator.py
2.) Put all lightcurve outputs from xray_generator.py into a single directory
    * Note: run_stats will calculate an energy distribution and flare frequency
      using *all* lightcurves in the specified directory. This allows the user 
      to generate a single energy distribution using any number of lightcurves
3.) Update the "Input Parameters" (labeled in the script) to match those
    used in xray_generator.py 
4.) Run the script!  

################################

Outputs: 
  - A plot showing the simulated energy distribution of flares, compared
    to a target or observed energy distribution
    - if plot_flare_freq = True then the flare frequency for each model
      will also be plotted, compared to the target flare frequency 

  - Chi**2: represents how well the modeled data matches the observed data
    defined as: 
         Chi**2 = (Simulated Distribution - Target Distribution)**2 / Target Distribution
    Smaller Chi**2 values indicate a better fit.

  -Average number of flares: determines the average number of flares that occur for 
   all light curves in the simulated time.  If only one light curve is run, then 
   this is the total number of flares.
#################################
'''
import numpy as N
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import pdb as pdb
import csv as csv
import glob as glob
from matplotlib import rc
from datetime import datetime

# set to true for nicer looking text
rc('text',usetex=False)

beg_time = datetime.now()

#####################
# Function that creates a histogram
#####################
def histSpec(x,y):
   """
   Plot a histogram style boxy-plots.
   """
   np = len(x)
   xp = N.zeros(np*2 + 2, dtype=N.float)
   yp = N.zeros(np*2 + 2, dtype=N.float)
   for bb in range(np):
       if bb > 0:
           yp[2*bb + 0] = y[bb-1]
           xp[2*bb + 0] = x[bb-1]+(x[bb]-x[bb-1])/2.0
           xp[2*bb + 1] = x[bb-1]+(x[bb]-x[bb-1])/2.0
       else:
           xp[2*bb + 0] = x[bb]-(x[bb+1]-x[bb])/2.0
           xp[2*bb + 1] = x[bb]-(x[bb+1]-x[bb])/2.0
           yp[2*bb + 0] = 0.0
       yp[2*bb + 1] = y[bb]
   xp[np*2+0] = x[np-1]+(x[np-1]-x[np-2])/2.0
   yp[np*2+0] = y[np-1]
   xp[np*2+1] = x[np-1]+(x[np-1]-x[np-2])/2.0
   yp[np*2+1] = 0.0
   return (xp,yp)

#################
# This function plots the target energy distribution
# of flares. Note that this distribution currently uses 
# an observed flare energy distribution presented in 
# Wolk et al. (2005), see Figure 9
# Describes solar mass PMS stars
################
def observed_dist(Emin,Emax,Estep):
    freq = []
    energy = []
    alpha = 1.66
    for i in N.arange(Emin,Emax,Estep):
        f = 1.1 * (10**i / 10**Emin)**(-(alpha-1.0))
        freq.append(f)
        energy.append(10**i)
    return freq,energy

################
# Function find_flares:
#     Search for individual, distinguishable flares from a light curve
#     Individual flares are constructed by:
#       1. Search for local maxima (find each flare)  
#          A point, P, is considered a maxima if 
#          a. that point is greater than P+1 and P-1
#          b. P/(P+1) and P/(P-1). Ensures P is a sharp peak (1.5% greater than P+1 and P-1)
#       2. find each flare beginning and end
#
#     Inputs:   
#           Lchar   : observed characteristic luminosity (erg/s) used in xray_generator.py
#           Emin    : minimum energy required to be classified as a flare (erg)
#           curve   : light curve produced by xray_generator.py
#                     units are total energy (erg) including Lchar
#           alltime : list containing every time step from the light curve
#     Outputs:  
#           Total flare energy : list containing the total energy (erg) of each distinguishable flare
#           totflares          : total number of distinguishable flares           
################
def find_flares(Lchar,Emin,curve,alltime,normalize):
    totflares = []
    Lchar = 10**Lchar
    min_flare = 0.1*Lchar # minimum flare energy considered
    min_end = 10**30.20   # minimum energy that defines end of a flare
    flare_energy   = [] # y axis

    timestep = alltime[1] - alltime[0] #assumes uniform timesteps

    curve = N.array(curve)
    alltime = N.array(alltime)
   
    # 1. Search for local maxima
    ratio_left = curve[1:] / curve[:-1]
    ratio_right = curve[:-1] / curve[1:]
    ratio_left = N.append([0],ratio_left)
    ratio_right = N.append(ratio_right,[0])

    min_crit = curve >= min_flare
    left_crit = ratio_left > 1.015
    right_crit = ratio_right > 1.015
    peaks_array = min_crit*left_crit*right_crit
    peaks = N.where(peaks_array)
    peaks = peaks[0]
    slopes = (curve[1:] - curve[:-1]) / (alltime[1:] - alltime[:-1])

    rat = 0.95 # minimum percent decrease from flare peak allowed to end the flare

    # 2. find the beginning and end of each flare
    for cnt in peaks:
        point = curve[cnt]
        
        time_after = alltime[cnt:] # time steps after flare peak
        time_before = alltime[:cnt] # time steps before flare peak

        # determine the end of the flare (decay)
        future_ratio = (curve[cnt:]/point)<rat #all future ratios
        future_low_lum = curve[cnt:] <= min_end #find all future points that above the min luminostiy 
        future_slopes = slopes[cnt-1:] > 0.0 #fining all future pointss that have a positive slope
        flare_stopper = N.zeros(len(time_after),dtype=bool)
        flare_stopper[-1] = True       

        end_flare_sum = future_low_lum + (future_ratio * future_slopes) +flare_stopper#if either flare1 or flare2 are true, then the sum is true
        end_flare = N.where(end_flare_sum) #an aray of indices
        end_flare = end_flare[0]
        eflare = curve[cnt:cnt+end_flare[0]] 
        etimes = alltime[cnt:cnt+end_flare[0]]

        # determine beginning of the flare (rise)
        past_ratio = (curve[:cnt]/point)<rat
        past_low_lum = curve[:cnt] <= min_end #10**30.24
        past_slopes = slopes[:cnt] < 0.0
        
        flare_stopper = N.zeros(len(time_before),dtype=bool)
        flare_stopper[0] = True
        beg_flare_sum = past_low_lum + (past_ratio * past_slopes) + flare_stopper
        beg_flare = N.where(beg_flare_sum)
        beg_flare = beg_flare[0]
        flare_beg = beg_flare[-1] +1
        bflares = curve[flare_beg:cnt]
        btimes = alltime[flare_beg:cnt]
        
        # construct the individual flare
        ind_times = N.append(btimes,etimes)
        ind_lum = N.append(bflares,eflare)
        ind_flare = ind_lum - Lchar
        energy = N.trapz(ind_flare,ind_times) #total energy, erg
        

        if energy >= 10**Emin: # only considered a flare if the total energy is greater than Emin
            if normalize == True: 
                ind_lum = N.array(ind_lum)/Lchar
            else: 
                pass
            plt.scatter(N.array(ind_times)/86400.0,ind_lum,marker='x')
            #plt.plot(N.array(ind_times)/86400.0,ind_lum)
            totflares.append(time[cnt])
            flare_energy.append(energy)
        else:
            pass
    return N.log10(flare_energy),totflares

#####################
# Function read_curve: 
#     Read in the lightcurve generated by xray_generatory.py
# Inputs: 
#     filename: name of the lightcurve
# Outputs: 
#     Lxr: Luminosity of the read in lightcurve at each time step
#     time: time steps for the light curve (seconds)
#####################
def read_curve(filename):
    days = 86400.0 # seconds in 1 day
    Lxr = []
    time = []

    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file,delimiter=",")
        line_cnt = 0
        for row in csv_reader:
            if row[0] == "#":
                pass
            else:
                time.append(float(row[0]))
                Lxr.append(float(row[1]))
    
    plt.plot(N.array(time)/days,Lxr,label=r'light curve',color='black')
    return Lxr,time

############################################
############################################
#
# INPUT PARAMETERS
#
############################################
############################################
# normalized_curve  if the lightcurve output was normalized to deltaL
#                   instead of a total energy output (Etot, erg),
#                   set equal to True. Else set to False
# show_curve        will plot each indiviudal curve showing each flare
#                   that was found to be observable and distinguishable
#                   each flare will be a different color 
# normalize_edist   when True, will normalize the energy distribution curve to 1, 
#                   so that 1 flare occured at the minimum observable energy
#                   when False, plots total number of flares of each energy   
# plot_flare_freq   when True, the observed flare frequency for each light curve
#                   will be plotted.
#                   note: is not useful when only one light curve is generated
###################
min_seed = 1
max_seed = 50
Lchar_obs = 10**30.25 #erg/s, observable characteristic luminosity
lchar = 10**30.25     # lchar used in xray_generatory.py to create the lightcurve
flare_freq = 49.0     # target average number of flares that occur in 1 year
Emin = 34.0           # log(Emin) erg, minimum total energy of an *observable* flare
                      # note that Emin_modeled is likely < Emin_observable
Emax = 37.57          # log(Emax) erg, maximum total energy allowed for a flare
Estep = 0.1           # Energy time step resolution 
alpha = 3.4           # powerlaw index, absolute value

normalized_curve = True
show_curve = False
normalize_edist = True

plot_flare_freq = True

model_dir = './'      # location of all lightcurve outputs 
target_data = './wolk05_edistdata.csv'
print('NOTE: default target energy distribution is taken from Wolk et al. (2005) Figure 9. This energy distribution is for a solar mass premain sequence star.')
##########################################
##########################################

lchar = N.log10(lchar)

filelist = glob.glob(model_dir+'lightcurve*') 

Edist = []
flr_cnt = []

for f1 in filelist:

    labl = 'Light curve'
    lightcurve,time = read_curve(f1)

    if normalized_curve == True:
        lightcurve = N.array(lightcurve) * Lchar_obs
    else:
        pass

    flare_energy,resolved_peaktime = find_flares(lchar,Emin,lightcurve,time,normalized_curve)
    #print("Total Number of Flares: %i"%(len(resolved_peaktime)))
    
    if show_curve == True:

        plt.legend()
        plt.xlabel(r"time (days)")
        if normalized_curve == True: 
            plt.ylabel(r"$\Delta L_{\rm XR}$")
        else: 
            plt.ylabel(r"E$_{\rm tot}$ (erg)")
        plt.show()
    
    else:
        plt.clf()
    
    Edist = Edist + flare_energy.tolist()
    flr_cnt.append(len(resolved_peaktime))

#
# Plot target (observed) energy distribution and simulated energy distribution
#
if plot_flare_freq == True:
    plt.subplot(121)
else: 
    plt.subplot(111)
plt.xlabel(r"E$_{\rm tot}$ (erg)")
if normalize_edist == True:
    plt.ylabel(r"Cumulative Number of Flares (normalized)")
    obs_freq,obs_en = observed_dist(Emin,Emax,Estep)
    plt.semilogx(obs_en,obs_freq,label=r"Target Energy Distribution, Wolk+(2005)",color='grey')
else:
    plt.ylabel(r"Cumulative Number of Flares (not normalized)")

numbins = 29
Evec = N.linspace(Emin,Emax,numbins)
num_per_bin,binval = N.histogram(Edist,bins=Evec)

for ind in range(len(num_per_bin)):
    num_per_bin[ind] = sum(num_per_bin[ind:])

float_per_bin = [N.float(x) for x in num_per_bin]
if normalize_edist == True:
    histEx,histEy = histSpec(10**(binval[0:-1]),float_per_bin/(N.max(float_per_bin)))
else:
    histEx,histEy = histSpec(10**(binval[0:-1]),float_per_bin)#/(N.max(float_per_bin)))
plt.plot(histEx,histEy,label=r'Modeled E$_{{\rm dist}}$',color='blue')

energy_data = [10**Emin]
freq_data = [1.0]
with open(target_data) as csv_file:
    csv_reader = csv.reader(csv_file,delimiter=',')
    for row in csv_reader:
        if row[0][0] == "#":
            pass
        else: 
            en = 10**(float(row[0]))
            num = float(row[1])
            energy_data.append(en)
            freq_data.append(num)
energy_data.append(10**Emax)
freq_data.append(1e-10)

#
# interpolate a curve for the simulated data from the model
#
float_per_bin = N.array(float_per_bin)/max(float_per_bin) #normalize 
evec = 10**Evec

f_data = interp.interp1d(energy_data,freq_data,fill_value='extrapolate')
f_sim  = interp.interp1d(evec[:-1],float_per_bin,fill_value='extrapolate')

allE = N.arange(Emin,Emax,Estep)
allE = 10**allE

freqdata = f_data(allE)
freqsim  = f_sim(allE)

chi = (freqdata-freqsim)**2/(freqdata+freqsim)
    
chi2 = sum(chi)
print("Chi**2 = %.5f"%(chi2))

plt.legend()

#
# Plot the Flare Frequency for each Seed
#
if plot_flare_freq == True:
    plt.subplot(122)

    seed_list = list(range(min_seed,max_seed))
    flist = [1.0]*(len(seed_list))
    norm_flare_freq = flare_freq*time[-1]/3.154e7
    freq_upper = 1.1 * norm_flare_freq * N.array(flist)
    freq_lower = 0.9 * norm_flare_freq * N.array(flist)

    plt.plot(seed_list,freq_upper,color='grey',linestyle='-')
    plt.plot(seed_list,freq_lower,color='grey',linestyle='-',label='+/- 10% of target flare frequency')
    plt.scatter(seed_list,flr_cnt,color='red',marker='x')

    plt.ylabel(r"Total Number of Observable Flares")
    plt.xlabel(r"Model Number")
    plt.legend()


end_time = datetime.now()
avg = sum(flr_cnt)/len(flr_cnt) #average number of flares
print("Average Number of Flares: %.2f per %.2f sec"%(avg,time[-1]))
print("Began: %s "%(beg_time))
print("Ended: %s "%(end_time))

#plt.tight_layout()

#plt.savefig('flare_stats.pdf',dpi=None,facecolor='w',edgecolor='w',orientation='portrait',format='pdf',transparent=False)
plt.show()


##################################
# End of run_stats.py
##################################
