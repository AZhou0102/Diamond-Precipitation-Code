#!/usr/bin/env python
# coding: utf-8

# In[107]:


#Calculate rho(r)n+1 from P(r)n using equation of state
#Units all in SI

import numpy as np

#computes pressure from density via murnaghan. Used for checking density guesses
def murnaghan(rho, rho_0, B_0, B_01):
    pressure = (B_0/B_01)*((rho_0/rho)**(-1*B_01)-1)
    return(pressure)

#takes a target pressure and returns the density associated with it via murnaghan
def invert_murnaghan(P_target, rho_0, rho_min, rho_max, B_0, B_01, allowed_discrepancy = (0.05 * (10 ** 9))):
    rho_guess = (rho_min + rho_max) / 2
    P_guess = murnaghan(rho_guess, rho_0, B_0, B_01)

#    print('rho guess: ', rho_guess)
#    print('guess: ', P_guess)
#    print('target: ', P_target, '\n')
    
    if abs(P_guess - P_target) < allowed_discrepancy:
#        print(f"rho_guess: {rho_guess} in kg/m^3,\nP_guess: {np.format_float_scientific(P_guess, precision = 9)},\nP_target: {np.format_float_scientific(P_target, precision = 9)} in Pa\n")
        return(rho_guess)

    elif P_guess > P_target:
        return(invert_murnaghan(P_target, rho_0, rho_min, rho_guess, B_0, B_01, allowed_discrepancy))

    elif P_guess < P_target:
        return(invert_murnaghan(P_target, rho_0, rho_guess, rho_max, B_0, B_01, allowed_discrepancy))


#example list of target pressures
example_target_pressures = [100000000000, 500000000000, 900000000000]
wolanin_97 = [1487, 1487, 14870, 14.9 * (10 ** 9), 6.2]

#turns list of pressures and list of parameters into an iterable dictionary to unpack for inputs
def make_input(target_pressures, parameters):
    inputs = {}
    for i in target_pressures:
        inputs[i] = [i] + parameters
    return(inputs)

#list parameters in order of: rho_0, rho_min, rho_max, B_0, B_01
def display_invert_murnaghan(target_pressures, parameters):
    for key, value in make_input(target_pressures, parameters).items():
        print(invert_murnaghan(*value))

#example
#display_invert_murnaghan(example_target_pressures, wolanin_97)

ambient = [0, 50.47 * (10 ** 9), 128* (10 ** 9)]
stichorite = [2143.5, 2143.5, 21435, 309.9 * (10 ** 9), 4.59]
#display_invert_murnaghan(ambient, stichorite)

