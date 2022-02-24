'''
#################################
X-ray Light Curve Generator: XGEN
################################# 

This model was written by Abygail R. Waggoner(1*) and L. Ilsedore Cleeves(1,2).

XGEN generates a light curve based on a set of parameters that can be determined
by observations. XGEN was inspired by the statistical analysis of X-ray flaring events described 
in Wolk et al. (2005) and the mathematical analysis of X-ray flaring events
described in Caramazza et al. (2007). 

Note that this model was originally written to produce an X-ray flare lightcurve, 
but it can be used to produce any stochastic light curve based on a power law

(1) University of Virginia, Department of Chemistry
(2) University of virginia, Department of Astronomy
* contact information: arw6qz@virginia.edu

###########
# OVERVIEW 
###########

This script generates a light curve based on the power-law 

dN/dE = beta * logE ** -alpha

where N is the number of flares that over the energy range dE, 
        E     = the total energy of the flare
        beta  = a normalization constant dependent on average flare frequency
        logE  = the total integrated energy of the flaring event
        alpha = the power law index

To run XGEN: go to the bottom of the script,
where "Variables/Physical Parameters" are defined

#########
# GENERATED LIGHT CURVES
#########

Lightcurves are generated using a random number generator, 
where the probability (P) of a flare within a set energy range,
   P = beta * (E2**(-alpha+1) - E1**(-alpha+1))
   beta = -1.0 * flare_frequency * time_step * Emin**(alpha)
occuring is determined by the power law above.
The total energy of each flare is set as a bin, where the upper bound for 
Etot is set by E2, and the lower bound is set by E1. The size of the 
energy bin is set by deltaE = E2 - E1

All flares are assumed to have an exponential rise and decay time,
set by trise and tdecay. Currently, trise and tdecay are uniform for all 
flares, but the model can be updated so that tirse and tdecay are determined
based on a distribution. 

The lightcurve is constructed by 
   1) determing at which time steps flares of energy Etot occur
      this is determined by the random number genorator/probability distribution
   2) The maximum luminosity (peak) of each flare is determined by 
      Etot = Lpeak * (flare_rise + flare_decay)
      flare_rise = integral(exp(t/trise)dt) from -inf to tpeak
      flare_decay = integral(exp(-t/tdecay)dt) from tpeak to inf
   3) A final, integrated light curve is constructed by summing over all flares
      at each individual time step

This model can generate a series of light curves with random seeds ranging from
seedmin to seedmax. The series of light curve outputs can be read into 
run_stats.py, which will plot the energy distribution of all *observable*
flares and the number of *observable* flares at each seed 
Observable flares are set by the user in run_stats.py.


##########
# OUTPUT FILES
##########
   
    light_curve_"file_name".csv      
       X-ray flare light curve at resolution tstep (seconds)
       file_name is defined within the output file
         and includes all parameters used to generate 
         the light curve
       The light curve can be either
         1: total energy at resolution tstep
            (energy of flares + baseline/characteristic energy)
            set by normalize = False
         2: relative change in energy when there is a flare, versus
            when there is not
            DeltaL = (L_flares + L_characteristic) / L_characteristic 
            set by normalize = True
        This file can be read into 
         1: plot_curve.py
            plots luminosity over time
         2: run_stats.py 
            identify *observable* flares and plot flare energy distribution
            and frequency 
           
    peaks_"file_name".inp                   
       includes all individual flare peaks (DeltaL) from the light curve, 
       time of peak (seconds), and flare rise and decay times (seconds)

    xgen_parameters.inp
        creates a file containing all of the input parameters used to 
        generate the light curve(s)

'''
import csv
import numpy as N
import pdb as pdb
from datetime import datetime
import time as comptime

starttime = datetime.now()

############################
# Function: create_xparams_output
#    Creates an output file contain all parameters
#    use to generate the light curve
# Inputs: 
#    all input parameters are defined under
#    "VARIABLES AND PHYSICAL PARAMETERS"	
# Outputs: 
#    a .inp file 
#    nothing returned
############################ 
def create_xparams_output(timelapse,timestep,seedmin,seedmax,flare_freq,Emaximum,\
                          Eminimum,Emin_observed,E_step,t_rise,t_decay,alpha,\
                          Lcharacteristic,Lchar_obs):
    file_name = 'xgen_parameters.inp'
    f0 = open(file_name,'w+')

    f0.write('# Physical parameters used in XGEN to produce the light curve \n')

    f0.write('timelapse = %.1f \n'%(timelapse))
    f0.write('timestep = %.2f \n'%(timestep))
    f0.write('seed_range = %i to %i\n'%(seedmin,seedmax-1)) #assumes only 1 iteration was ran
    f0.write('avg_flr_freq = %.6f \n'%(flare_freq))
    f0.write('Emax = %.2f \n'%(Emaximum))
    f0.write('Emin = %.2f \n'%(Eminimum))
    f0.write('Emin_obs = %.2f \n'%(Emin_observed))
    f0.write('Estep = %.3f \n'%(E_step))
    f0.write('# trise and tdecay are currently uniform for all flares\n')
    f0.write('trise = %.2f \n'%(t_rise))
    f0.write('tdecay = %.2f \n'%(t_decay))
    f0.write('# dN/dE = kE**-alpha\n')
    f0.write('alpha = %.2f \n'%(alpha))
    f0.write('Lchar = %.2f \n'%(Lcharacteristic))
    f0.write('Lchar_obs = %2.f \n'%(Lchar_obs))
    f0.close()

    return

############################
# Function create_Lmax_output:
#    Creates an output file including the peak luminosity of
#    each flare greater than the observed characteristic luminosity
#    This function is set to always output the normalized flare peaks, 
#    but that can be changed by setting "normalize = False" inside 
#    this function
#    This output is meant to be read into the chemical disk model
# Inputs:
#    trise:      exp. rise time of flares (sec) (currently uniform for all flares)
#    tdecay:     exp. decay time of flares (sec) (currently uniform for all flares)
#    Lchar:      characteristic luminosity (erg/s) used in calculations. Does not include microflares
#    Lchar_obs:  observed characteristic luminosity (erg/s). Includes all microflares
#    Lpeak:      list of the peak luminosity of each flare (erg/s)
#    t_max:      list of the time each flare peak occurs (sec)
#    filename:   name of the output file. Contains physical parameters of the model
#    name_order: a string defining the order of variables in filename
# Output:
#    a .inp file containing a tpeak, Lpeak, trise, and tdecay
#    nothing retured
############################
def create_Lmax_output(trise,tdecay,Lchar,Lchar_obs,Lpeak,t_max,filename,name_order):    
  
    if normalize == True:
        Lpeak = (N.array(Lpeak)+Lchar)/Lchar
        out_type = 'LXR is xray energy compared to characterstic xrays -> LXR = (Eflare + Lchar)/Lchar'
    else:
        Lpeak = N.array(Lpeak) + Lchar
        out_type = 'LXR is the total energy (Lflare + Lcharacterisctic) of a flare peak over time (erg/s)'
   
    spc = "     " # Evenly spaced columns in output file
    
    # create file with every flare peak
    f = open("peaks_%s.inp"%(filename),"w+")
    
    f.write("# Produced by xray_generator.py \n")
    f.write("# Parameters: %s = %s\n"%(name_order,filename))
    f.write("# tpeak (s)           Lpeak (erg/s)          Trise (s)  Tdecay (s)\n")
    
    for peaktime in t_max:     
        lumpeak = Lpeak[t_max.index(peaktime)]
        f.write("%.9e%s%.9e%s%.2f%s%.2f\n"%(peaktime,spc,lumpeak,spc,trise,spc,tdecay)) 
    
    f.close()

    return 

#########################
# Function create_curve_output:
#    Creates an output file of the integrated light curve over time
#    at a resolution set by timestep.
#    This ouput is meant to be read into plot_curve.py and run_stats.py
# Inputs:
#    time:       list of every time step (sec)
#    Lchar:      characteristic luminosity (erg/s) used in calculations. Does not include microflares
#    Lchar_obs:  observed characteristic luminosity (erg/s). Includes all microflares
#    light_curve:list of X-ray luminosity at every time step (erg/s)
#    filename:   name of the output file. Contains physical parameters of the model
#    name_order: a string defining the order of variables in filename
#    normalize:  boolian that determines if the output luminosity is in total energy or 
#                normalized energy (DeltaLXR)
# Output:
#    a .csv file containing light_curve and time
#    nothing retured
########################    
def create_curve_output(time,Lchar,Lchar_obs,light_curve,name_order,flare_count,normalize,filename):

    # output light curve at resolution timestep    
    if normalize == True: #output is energy with respect to characteristic luminosity (unitless)
        cnt = 0
        for i in light_curve: 
            light_curve[cnt] = (i+Lchar)/(10**Lchar_obs)
            cnt += 1
        out_type = 'Lightcurve is xray strength compared to\
                     characteristic X-rays -> DeltaLXR = (Lflare + Lchar)/Lchar'
    else: # output is total energy (ergs), including characteristic
        cnt = 0
        for i in light_curve:
            light_curve[cnt] = i+Lchar
            cnt += 1
        out_type = 'Lightcurve is total energy (ergs) of X-ray flares plus the characteristic energy'
    
    f = open("lightcurve_%s.csv"%(filename),"w+")
    
    # make sure there is a comma after #, 
    #otherwise the plotting script won't be able to distinguish the comments 
    f.write("#, Light curve produced by xray_generator.py\n")
    f.write("#, File parameters are: %s = %s\n"%(name_order,filename))
    f.write("#, Columns are comma separated\n")
    f.write("#, time(s),LightCurve\n")

    for t in time: 
        l = light_curve[time.index(t)]
        f.write("%f,%f\n"%(t,l)) 
    f.close()

    return

##########################################
# Function generate_light_curve:
#    From a list of the total energy of each flare (produced by generate_xrays), 
#    first find the peak luminosity of each flare, 
#    then construct a single integrated light curve at resolution timestep.
#    Light curve construction is described at the top of this script/in the user manual
# Inputs:
#    time:       list of every time step (sec)
#    t_max:      list of the time each flare peak occurs (sec) (produced in xray_generator)
#    Lchar:      characteristic luminosity (erg/s) used in calculations. Does not include microflares
#    Etot:       list of the total energy of each individual flare (produced in xray_generator) 
#    trise:      exp. rise time of each flare (sec) (currently uniform for all flares)
#    tdecay:     exp. decay time of each flare (sec) (currently uniform for all flares)
# Returns:
#    Lpeak:      list of the peak luminosity of each individual flare (erg/s)      
#    light_curve: list of integrated X-ray luminosity at every time step (erg/s)     
##########################################
def generate_light_curve(time,t_max,Lchar,Etot,trise,tdecay):
    Lpeak = []
    t_max = N.array(t_max)
    time = N.array(time)
    time_step = time[1] - time[0]
    
    light_curve = [0] * len(time) 
    
    # define the beginning and end of each flare   
    # increase beg and end if you find the flare is being cut off to early
    # Note: increasing beg and end will slow down the code
    beg = 6.0 # num of exp. rise times considered before flare peak
    end = 25.0 # num of exp. decay times considered after flare peak 
    flare_beg = t_max - (beg*trise)
    flare_end = t_max + (end*tdecay)   
    
    e_count = 0
    for E in Etot: #for each flare
 
        # determine time steps and energy at each time step of the flare rise and decay
        rise_time = N.arange(flare_beg[e_count],t_max[e_count],time_step)
        prime_rise_time = rise_time - t_max[e_count]  #normalizes t_max to be zero 
        rise_energy = N.exp(prime_rise_time/trise)        
        decay_time = N.arange(t_max[e_count]+time_step,flare_end[e_count]+time_step,time_step)
        prime_decay_time = decay_time -t_max[e_count] #normalize t_max to be zero
        decay_energy = N.exp(-1.0*prime_decay_time/tdecay)
        
        # integrate over the flare rise and decay
        int_decay = N.trapz(decay_energy,prime_decay_time)
        int_rise = N.trapz(rise_energy,prime_rise_time)
    
        L_peak = (10**E)/(int_decay + int_rise)
        Lpeak.append(L_peak)        
 
        # determine the individual light curve at every time step and 
        # include the individual flare in the integrated light curve
        prime_time = time - t_max[e_count]
        peak_location = (prime_time.tolist()).index(0.0)
        flare_rise  = L_peak*(N.exp(prime_time[0:peak_location]/trise))
        flare_decay = L_peak*(N.exp(-1.0*prime_time[peak_location:]/tdecay))
        light_curve = light_curve + N.concatenate([flare_rise,flare_decay])
        
        e_count += 1

        if (e_count%100) == 0: 
            print("Curve %i of %i complete"%(e_count,len(Etot)))
        else: 
            pass

    return Lpeak,light_curve


##################################
# Produce Stochastic Flares with Total Energy E
# Function xrays_random:
#    Use a random number generator to determine whether or not a flare
#    between a total flare energy of E1 and E2 occurs at every time step
#    This is the main function of XGEN and calls all other functions
# Inputs: 
#    all input variables are defined in "Variable/Physical Parameters"
# Outputs: 
#    None. It calls all other functions
##################################
def xrays_random(index_alpha,randseed,time_lapse,time_step,flare_freq,Emax,\
                 Emin,Emin_obs,Estep,Lchar,Lchar_obs,trise,tdecay,filename,\
                 name_order,normalize):
    
    trise = trise * 3600.0 #converts to seconds
    tdecay = tdecay* 3600.0
    
    N.random.seed(randseed)
    
    Etot = []
    time = []
    t_max = []
    flare_prob = []
    energy_range = []

    ia = index_alpha -1.0

    # normalization constant
    # NOTE beta is currently divided by 2 to achieve the correct observed flare frequency
    beta = 0.5 * -1.0 * (flare_freq * time_step * (10.0**Emin_obs)**ia)

    # create an array of all possible total flare energies at resolution Estep
    energy_range = N.arange(Emax,Emin-Estep,-1.0*Estep)
    E1 = 10**(energy_range+Estep) 
    E2 = 10**(energy_range-Estep)
    flare_prob = beta * (E1**-ia - E2**-ia)
    NE = len(energy_range.tolist())

    time = N.arange(0.0,time_lapse+time_step,time_step) #array of all time steps
    NT = len(time)
    
    # create an array of the probability of each flare with E2 < Etot < E1 
    # at every time step
    flare_prob_arr, time_mesh = N.meshgrid(flare_prob,time)
    
    # identify every position in flare_prob_arr where a flare occured
    # a flare occurs everywhere the random numer is less than the flare probability 
    rns = N.random.sample([NT,NE])
    inds = N.where((flare_prob_arr - rns) > 0) 
    
    Etot = energy_range[inds[1]].tolist() #an array of all flares of energy Etot that occured
    t_max = time[inds[0]].tolist() #an array of the time of each flare peak
    time = time.tolist()

    print('random flares determined')
    
    # Generate the light curve
    Lmax,light_curve = generate_light_curve(time,t_max,Lchar,Etot,trise,tdecay)
    
    # Creates output files
    create_curve_output(time,Lchar,Lchar_obs,light_curve,name_order,len(Lmax),normalize,filename)
    create_Lmax_output(trise,tdecay,Lchar,Lchar_obs,Lmax,t_max,filename,name_order)

    return

##########################################
##########################################
#
# VARIABLES AND PHYSICAL PARAMETERS
# 
##########################################
##########################################

# Time
timelapse       = 31450000.0          # Simulation length, seconds
timestep        = 3600.0              # time resolution, in seconds 

# Seed of the random number generator
# if only 1 lightcurve is to be generated, set seedmax = seedmin + 1
seedmin         = 1
seedmax         = 50

# Physical Flare Parameters
flare_freq      = (1.0/650.0)/1000.0  # Observed number of flares per second
Emaximum        = 37.57               # Maximum log(total energy) of a flare (logE erg)
Eminimum        = 32.5                # Minimum log(total energy) of flares that can occur in the model (logE erg)   
Emin_observed   = 34.0                # Minimum log(total energy) of flares that are cosidered detectable
                                      # NOTE. Emin_observed changes the flare energy distribution and 
                                      # flare frequency (see beta in xrays_random()) 

E_step          = 0.01                # Energy resolution. Keep at/below 0.01
t_rise          = 3.0                 # Exponential rise time of the flare (hours)
t_decay         = 8.0                 # Exponential decay time of the flare (hours)
alpha           = 1.64                # absolute value, dN/dE = kE**-alpha

# Normalize the final ouput file with respect to the characteristic luminosity? 
# normalize = True: Output is deltaLxr = (Eflare + Lchar)/Lchar
# normalize = False: Output is total energy (ergs) Lxr = Lflare + Lcharacteristic
normalize = True

Lcharacteristic = 30.25               # Baseline/characteristc luminosity of the source (logL erg/sec)
                                      #    Note: May be less than the Lchar_obs
                                      #          Does not include micro or nanoflares
Lchar_obs = 30.25                     # average observed characterisitc luminosity (logL erg/sec)
                                      #    Note: Contains micro and nano flares
                                      # If uncertain, set Lcharacteristic = Lchar_obs 

#########################################
#########################################
#########################################

observed_flarefreq = timelapse*flare_freq
seed_vals = []
tot_flare_cnt = []
obs_flare_cnt = []
E_stats = []

for i in range(seedmin,seedmax,1):

    rand_seed = i

    #if file_name is changed, also update name_order
    file_name = "%.2f_%.1f_%.1f_%.2f_%.2f_%.2f_%.2f_%.2f_%.1f_%.1f_%.2f_%i"%\
                (timelapse/31556952.0,timestep/60.0,flare_freq*1000.0,\
                 Emaximum,Eminimum,E_step,Lcharacteristic,Lchar_obs,t_rise,t_decay,alpha,rand_seed)
    name_order = "simlength_timestep_flarefreq_Emax_Emin_Estep_Lchar_Lcharobs_trise_tdecay_alpha,seed"

    xrays_random(alpha,rand_seed,timelapse,timestep,\
                     flare_freq,Emaximum,Eminimum,Emin_observed,E_step,10**Lcharacteristic,\
                     Lchar_obs,t_rise,t_decay,file_name,name_order,normalize)
    
    print("Iteration #%i Completed"%(rand_seed))


create_xparams_output(timelapse,timestep,seedmin,seedmax,flare_freq,Emaximum,\
                      Eminimum,Emin_observed,E_step,t_rise,t_decay,alpha,\
                      Lcharacteristic,Lchar_obs)

endtime = datetime.now()
timediff = endtime - starttime

print('Start time: %s'%(starttime))
print('End time: %s'%(endtime))
print('Time it took to run the model: %s'%(timediff))

##############################
# End of XGEN
##############################
