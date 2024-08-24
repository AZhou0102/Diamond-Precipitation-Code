#!/usr/bin/env python
# coding: utf-8

# In[198]:


#All units in SI

import numpy as np

#computes pressure from density via bm3. Used for checking density guesses
def bm3(rho, rho_0, B_0, B_01):
    pressure = (3 / 2) * B_0 * (((rho / rho_0) ** (7/3)) - ((rho / rho_0) ** (5/3))) * (1 + (3/4 * (B_01 - 4)) * (((rho / rho_0) ** (2/3)) - 1))
    return(pressure)



def invert_bm3(P_target = 999999999999, rho_0 = 1487, rho_min = 1487, rho_max = 14870, B_0 = 14.9 * (10 ** 9), B_01 = 6.2, allowed_discrepancy = (0.05 * (10 ** 5))):
    rho_guess = (rho_min + rho_max) / 2
    P_guess = bm3(rho_guess, rho_0, B_0, B_01)

#    print('rho guess: ', rho_guess)
#    print('guess: ', P_guess)
#    print('target: ', P_target, '\n')
    
    if abs(P_guess - P_target) < allowed_discrepancy:
#        print(f"rho_calculated: {rho_guess} in kg/m^3\nP_guess: {np.format_float_scientific(P_guess, precision = 9)}\nP_target: {np.format_float_scientific(P_target, precision = 9)} in Pa\n")
        return(rho_guess)

    elif P_guess > P_target:
        return(invert_bm3(P_target, rho_0, rho_min, rho_guess, B_0, B_01, allowed_discrepancy))

    elif P_guess < P_target:
        return(invert_bm3(P_target, rho_0, rho_guess, rho_max, B_0, B_01, allowed_discrepancy))

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
def display_invert_bm3(target_pressures, parameters):
    print("target pressures:", target_pressures, '\n')
    for key, value in make_input(target_pressures, parameters).items():
        print(invert_bm3(*value))

#example
#display_invert_bm3(example_target_pressures, wolanin_97)

ambient = [0, 50.47 * (10 ** 9), 128* (10 ** 9)]
stichorite = [2143.5, 2143.5, 21435, 309.9 * (10 ** 9), 4.59]
#display_invert_bm3(ambient, stichorite)

P_example = [0, 100000000000, 500000000000, 900000000000]
iron = [8269, 8269, 82690, 163.4 * (10 ** 9), 5.38]
#display_invert_bm3(P_example, iron)

