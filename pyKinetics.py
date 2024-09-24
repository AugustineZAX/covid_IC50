#! /usr/bin/env python
#!!!!! This code is for the publication level

#-------------------------------------------------------------------------------
# . File      : pyKinetics.py (python3.9.13)
# . Program   : MolarisTools
# . Copyright : USC, Aoxuan(Augustine) Zhang (Aug 4, 2024)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------

help_message = """
This python script helps to do a numerical integration of a set of coupled ODEs to determine
the time evolution of different species in the solvent, it is tested for the followint kinetic model
E + I <-> EI1 <-> EI2
in which '<->' denotes a reversible process

The IC50 example code is also provided, feel free to modify it to your kinetics!
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

#------------------------------
#
# We first get Kinetics for the natrual substrate as a reference for the following IC50 calculations
#
#------------------------------

ks_fw  = 1e-2
ks_bw  = 0.143
ks_cat = 0.43

# Define the ODEs
# ODEs for preliminary enzyme test, the result will be the reference for IC50 calculation
def odes00(x,t):
    
    # kinetic data
#     ks_fw  = 0.2
#     ks_bw  = 11.3
#     ks_cat = 0.43
    
    # assign each ODE to a vector element
    E  = x[0]  # apo enzyme, [nM] so as the following elements
    S  = x[1]
    ES = x[2]
    P  = x[3]
    
    # define each ODE
    dEdt  = -ks_fw*E*S + (ks_bw + ks_cat)*ES
    dSdt  = -ks_fw*E*S + ks_bw*ES
    dESdt = ks_fw*E*S - (ks_bw + ks_cat)*ES
    dPdt  = ks_cat*ES
    
    return [dEdt, dSdt, dESdt, dPdt]

# Give initial conditions
x00 = [0.5, 20, 0, 0]

# Declare a time vector
t_end = 2000
step_size = 0.05
tot_step = int(t_end/step_size)
t = np.linspace(0, t_end, tot_step)
x = odeint(odes00, x00, t)
E          = x[:,0]
S          = x[:,1]
ES   = x[:,2]
P      = x[:,3]

# plot
plt.plot(t, E, color = 'red')
plt.plot(t, S, color = 'orange')
plt.plot(t, P, color = 'green')
plt.plot(t, ES, color = 'blue')
plt.xlim(0,50)
plt.ylim(0,0.5)
plt.xlabel('time[s]')
plt.ylabel('concentration[\u03bcM]')
print (P[int(50/step_size)])
vP0 = P[int(50/step_size)]
ratio = (P[int(50/step_size)]/50)/(P[int(20/step_size)]/20)
print (ratio)

#---------------------------------------------------
#
#    ODE Define Block for Incubation (odes01)
#            and IC50 simulation (odes02)
#
#---------------------------------------------------


# For convenience, all rate constants are assigned out of function definition block
# But be sure to modify it if other purposes are served

k1_fw = 2e-2
k1_bw = 1.5e-2
k2_fw = 3.31
k2_bw = 1.97
# ks_fw  = 0.2
# ks_bw  = 11.3
# ks_cat = 0.43

# Define the ODEs
# ODEs for preincubation
def odes01(x,t):
    
    # kinetic data
#     k1_fw = 1e-3
#     k1_bw = 0.034195405
#     k2_fw = 74
#     k2_bw = 7e-4
    
    # assign each ODE to a vector element
    E        = x[0]  # apo enzyme, [nM] so as the following elements
    I        = x[1]  # unbound ligand
    EInoncov = x[2]  # complex for the 1st reaction
    EIcov    = x[3]  # complex for the 2nd reaction
    
    # define each ODE
    dEdt        = -k1_fw*E*I + k1_bw*EInoncov
    dIdt        = -k1_fw*E*I + k1_bw*EInoncov
    dEInoncovdt = k1_fw*E*I - (k1_bw + k2_fw)*EInoncov + k2_bw*EIcov
    dEIcovdt    = k2_fw*EInoncov - k2_bw*EIcov
    
    return [dEdt, dIdt, dEInoncovdt, dEIcovdt]

# Define the ODEs
# ODEs for IC50 simulation at time t
def odes02(x,t):
    
    # kinetic data
#     k1_fw = 1e-3
#     k1_bw = 0.034195405
#     k2_fw = 74
#     k2_bw = 7e-4
#     ks_fw  = 1e-2
#     ks_bw  = 0.143
#     ks_cat = 0.43
    
    # assign each ODE to a vector element
    E        = x[0]  # apo enzyme, [nM] so as the following elements
    I        = x[1]  # unbound ligand
    EInoncov = x[2]  # complex for the 1st reaction
    EIcov    = x[3]  # complex for the 2nd reaction
    S        = x[4]
    ES       = x[5]
    P        = x[6]
    
    # define each ODE
    dEdt        = - k1_fw*E*I + k1_bw*EInoncov - ks_fw*E*S + ks_bw*ES + ks_cat*ES
    dIdt        = - k1_fw*E*I + k1_bw*EInoncov
    dEInoncovdt = k1_fw*E*I - (k1_bw + k2_fw)*EInoncov + k2_bw*EIcov
    dEIcovdt    = k2_fw*EInoncov - k2_bw*EIcov
    dSdt        = - ks_fw*E*S + ks_bw*ES
    dESdt       = ks_fw*E*S - (ks_bw + ks_cat)*ES
    dPdt        = ks_cat*ES
    
    return [dEdt, dIdt, dEInoncovdt, dEIcovdt, dSdt, dESdt, dPdt]


    
# Now let's automize the process and get I0 vs inhibition plot
# Define a list for initial I0 you want to calculate:
I0_start = -3.0         # log space is used, the real start will be base**(I0_start)
I0_end   = 5
I0_list = np.logspace(I0_start, I0_end, 50)

# Define a list of incubation time
t1_list = [30, 15*60, 3*3600]
color_list = ['darkmagenta', 'darkcyan', 'brown']
marker_list = ['o', 'v', 's']
legend_list = ['30s', '15min', '3h']

for i in range(0, len(t1_list)):
    
    # We initialize the inhibition list
    inhibition_list = []
    
    for j in range (0, len(I0_list)):
        
        #--------------------------------------------
        #
        #     This is the incubation part
        #
        #--------------------------------------------
        
        # Give initial conditions
        x01 = [0.5, I0_list[j], 0, 0]

        # Declare a time vector
        t_end01 = t1_list[i]
        step_size01 = 0.05
        tot_step01 = int(t_end01/step_size01)
        t01 = np.linspace(0, t_end01, tot_step01)
        x01 = odeint(odes01, x01, t01)
        E          = x01[:,0]
        I          = x01[:,1]
        EInoncov   = x01[:,2]
        EIcov      = x01[:,3]
        
        if t1_list[i] == 0:
            
            Et1 = 0.5
            It1 = I0_list[j]
            EInoncovt1 = 0
            EIcovt1 = 0
        else:
                
            Et1 = E[-1]
            It1 = I[-1]
            EInoncovt1 = EInoncov[-1]
            EIcovt1 = EIcov[-1]
        
        #--------------------------------------------
        #
        #     This is the IC50 part
        #
        #--------------------------------------------
        
        # Give initial conditions
        x02 = [Et1, It1, EInoncovt1, EIcovt1, 20, 0, 0]

        # Declare a time vector
        t_end02 = 51
        step_size02 = 0.05
        tot_step02 = int(t_end02/step_size02)
        t02 = np.linspace(0, t_end02, tot_step02)
        x02 = odeint(odes02, x02, t02)
        E          = x02[:,0]
        I          = x02[:,1]
        EInoncov   = x02[:,2]
        EIcov      = x02[:,3]
        S          = x02[:,4]
        ES         = x02[:,5]
        P          = x02[:,6]
        
        inhibition = 100*(1 - P[int(50/step_size02)]/vP0)
        inhibition_list.append(inhibition)
#         print ('inhibition = ', round(inhibition, 2))
    
    if len(I0_list) == len(inhibition_list):
        plt.xscale ('log')
        plt.plot (1e3*I0_list, inhibition_list, color=color_list[i], marker=marker_list[i], 
                  label=legend_list[i], markersize=5)
        plt.ylim(0, 100)
#         plt.title('Nitrile Direct')
        plt.xlabel('[I] (nM)', font='arial', size=14, weight='bold')
        plt.ylabel('% Inhibition', font='arial', size=14, weight='bold')
        plt.legend()
        plt.axhline (y=50, color='grey', linestyle='dashed')
        
        plt.figtext(0.15, 0.77, '(a) Nitrile Direct', font='arial', size=12, weight='bold')

        
# plt.plot(E)
# plt.xlim(0,100)

#         plt.axvline (x=0.34*1e3, color='darkmagenta', linestyle='dashed')
#         plt.axvline (x=0.56*1e3, color='darkcyan', linestyle='dashed')
#         plt.axvline (x=0.76*1e3, color='brown', linestyle='dashed')
        
# Save the plot to the specified path with the given filename
path = '/Users/augustinedss/Desktop/0__DrWarshel/ashim_covid/matplotlib_sep05'
filename = 'Nitrile Direct'
full_path = f'{path}/{filename}.png'  # You can change the extension to .pdf, .svg, etc.

plt.savefig(full_path, dpi=1000)

#----------------------
#
#    For Automatic IC50 Searching with Given Kinetics
#
#----------------------

# Define Initial Conditions
I0 = 1e8                                 # Initial Guess of IC50, can be checked with a pre-run
time_list = [30, 15*60, 3*3600]         # time to measure IC50

for i in range(0, len(time_list)):
    
    inhibition = 100
    converge_check = abs(round(inhibition,2) - 50)
    I0_upper = I0
    
    
    while converge_check >= 0.1:

        if inhibition - 50 > 0:
            I_upper = round(I0, 3)
            I0 = round(I0/2, 3)
        elif inhibition - 50 < 0:
            I0 = round(0.5*(I_upper + I0), 3)

#         print ('I0 = ', I0)
        
        # Give initial conditions
        x01 = [0.5, I0, 0, 0]

        # Declare a time vector
        t_end = time_list[i]
        step_size = 0.05
        tot_step = int(t_end/step_size)
        t = np.linspace(0, t_end, tot_step)
        x = odeint(odes01, x01, t)
        E          = x[:,0]
        I          = x[:,1]
        EInoncov   = x[:,2]
        EIcov      = x[:,3]

        # # plot
        # # plt.plot(t, E, color = 'red')
        # # plt.plot(t, I, color = 'orange')
        # plt.plot(t, EInoncov, color = 'orange')
        # plt.plot(t, EIcov, color = 'green')
        # # plt.ylim(0,1e-2)
        # plt.xlabel('time[s]')
        # plt.ylabel('concentration[\u03bcM]')
        # print (E[-1], I[-1], EInoncov[-1], EIcov[-1])
        Et1 = E[-1]
        It1 = I[-1]
        EInoncovt1 = EInoncov[-1]
        EIcovt1 = EIcov[-1]

        # Give initial conditions
        x02 = [Et1, It1, EInoncovt1, EIcovt1, 20, 0, 0]

        # Declare a time vector
        t_end = 55
        step_size = 0.05
        tot_step = int(t_end/step_size)
        t = np.linspace(0, t_end, tot_step)
        x = odeint(odes02, x02, t)
        E          = x[:,0]
        I          = x[:,1]
        EInoncov   = x[:,2]
        EIcov      = x[:,3]
        S          = x[:,4]
        ES         = x[:,5]
        P          = x[:,6]
        Etot = E + EInoncov + EIcov +ES
        inhibition = 100*(1 - P[int(50/step_size)]/vP0)
        converge_check = abs(round(inhibition,2) - 50)

    print ('IC50(', time_list[i], ') = ', I0, '\n',
          'with [ES(', time_list[i], ')] = ', round(ES[-1],3), 'μM\n',
          'and inhibition = ', round(inhibition,2), '\n')

        # plot
        # plt.plot(t, E, color = 'red')
        # # plt.plot(t, I, color = 'orange')
        # plt.plot(t, EInoncov, color = 'green')
        # plt.plot(t, EIcov, color = 'turquoise')
        # # plt.plot(t, S, color = 'dodgerblue')
        # plt.plot(t, ES, color = 'darkviolet')
        # plt.plot(t, P, color = 'magenta')
        # # plt.plot(t, Etot, color = 'black')
        # # plt.xlim(0,300)
        # plt.ylim(0.45,0.5)
        # plt.xlabel('time[s]')
        # plt.ylabel('concentration[\u03bcM]')
        # print (E[-1], I[-1], EInoncov[-1], EIcov[-1], S[-1], ES[-1], P[-1])
        # print (P[int(50/step_size)])
    #         print ('inhibition = ', round(100*(1 - P[int(50/step_size)]/vP0),2),'%')
