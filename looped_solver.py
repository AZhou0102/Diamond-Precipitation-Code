#!/usr/bin/env python
# coding: utf-8
"""
Planetary Structure Solver

This module calculates density and pressure profiles inside a planet given a list of radii and
a dictionary describing its layered composition.

Inputs:
- radii_list: A list of radii (in meters) to evaluate.
  Recommended form: [*range(1, int(r_planet), int(r_planet / slices))]

- insert_dict: Dictionary from planetary_dictionary(planet_size, iron_part, rock_part)
  where:
    • planet_size is in Earth radii
    • iron_part is the fractional depth of the iron core
    • rock_part is the fractional depth of the rock mantle

Note:
- The largest radius in radii_list must match planet_size (converted to meters).
- All computations use SI units.

Returns:
- Two lists: densities and pressures corresponding to each radius.
"""

# In[6]:


import EoS_Bits as EOS
import planetary_dictionary as dct
import solve_adams_williamson as aw

import numpy as np
from matplotlib import pyplot as plt

r_earth = 6370 * (10 ** 3)


def Solver(radii_list, insert_dict, density_list = [], discrepancy = 10, calls = 0):
    '''
    Computes the self-consistent radial density and pressure profiles of a planet given a list of radii 
    and a planetary composition dictionary.

    Parameters:
    radii_list (list): A list of radial distances from the center of the planet (in meters) at which density and pressure will be calculated.
    insert_dict (dict): A dictionary specifying the planet's compositional layers. Can be created using planetary_dictionary() or manually.
    density_list (list, optional): An optional list of initial density values as guesses (eg: zero-pressure densities). Defaults to an empty list.
    discrepancy (float, optional): The convergence threshold for density in kg/m^3. Defaults to 10.
    calls (int, optional): Recursion counter to prevent infinite loops. Defaults to 0.

    Returns:
        tuple: A tuple containing two lists:
            - List of densities at each radius (in kg/m^3)
            - List of pressures at each radius (in Pa)
    '''
#The density list is initially empty and is populated after the first run.
#The discrepancy represents the error in density.


    pressure_n = {}
    pressure_n1 = {}
    density_n1 = {}
    calls += 1


    while calls <= 2000:
    #Prevents the number of recursions from diverging.
        if int(max([*insert_dict.keys()])) != int(max(radii_list)):            
            insert_dict[int(float(f"{int(max(radii_list))}"))] = insert_dict.pop(int(float(f'{max([*insert_dict.keys()])}')))
            return(Solver(radii_list, insert_dict, density_list = [], discrepancy = discrepancy, calls = 0))
        #Ensures that the maximum radius in the planet's dictionary and the radius list are compatible

        
        if len(density_list) == 0:
            for count in range(len(radii_list)):
                for cutoff in [*insert_dict.keys()]:
                        if radii_list[count] <= cutoff:
                            if count == len(density_list):
                                density_list.append(int(insert_dict[cutoff][0]))
        #On the first iteration, the Solver makes a density_list using initial density guesses as preliminary values
    
        pressure_list = aw.adams_williamson(radii_list, density_list)
        #pressure_list is a list of initial pressures at each radius.
        
        for n in radii_list:
            pressure_n[n] = pressure_list[n][1]
        #At each radius in the list, a corresponding pressure is assigned. This is the nth pressure.
        
        for n in radii_list:
            nth_p_n = pressure_n[n]
            #This assigns the pressure at the nth radius calculated for the nth time.
            for cutoff in list(insert_dict.keys()):    
                if n <= cutoff: 
                #Here, we select the correct EoS and EoS parameters corresponding that layer's radius, pressure, and material.
                    
                    if list(insert_dict.keys()).index(cutoff) == 0:
                    #This "if block" is for calculations within the Core.
                        
                        #We calculate the n+1th density and pressure using functions from EoS_Bits 
                    
                        density_n1[n] = EOS.CoreDensity(nth_p_n)
                        pressure_n1[n] = EOS.Vinet(density_n1[n],7678,136.2 * (10 ** 9),5.97)
                        break
                        
                    elif list(insert_dict.keys()).index(cutoff) == 1:
                    #This "if block" is for calculations within the Mantle.
                        
                        #We calculate the n+1th density and pressure using functions from EoS_Bits 
                        
                        density_n1[n] = EOS.RockDensity(nth_p_n)
                        
                        #Depending on the pressure, different EOS parameters are selected.
                        if nth_p_n < 2.5 * (10 ** 9):
                            pressure_n1[n] = 0.701 * EOS.BM3(density_n1[n], 3221,125 * (10 ** 9),4) + 0.299 * EOS.BM3(density_n1[n],2648,37.4 * (10 ** 9),6.2)
                            break
                        elif nth_p_n < 8 * (10 ** 9):
                            pressure_n1[n] = 0.701 * EOS.BM3(density_n1[n], 3221,125 * (10 ** 9),4) + 0.299 * EOS.BM3(density_n1[n],2921,96 * (10 ** 9),8.4)
                            break
                        elif nth_p_n < 14 * (10 ** 9):
                            pressure_n1[n] = 0.701 * EOS.BM3(density_n1[n], 3221,125 * (10 ** 9),4) + 0.299 * EOS.BM3(density_n1[n],4290,309.9 * (10 ** 9),4.59)
                            break
                        elif nth_p_n < 18 * (10 ** 9):
                            pressure_n1[n] = 0.701 * EOS.BM3(density_n1[n],3491,160 * (10 ** 9),4) + 0.299 * EOS.BM3(density_n1[n],4290,309.9 * (10 ** 9),4.59)
                            break
                        elif nth_p_n < 23 * (10 ** 9):
                            pressure_n1[n] = 0.701 * EOS.BM3(density_n1[n],3548,182 * (10 ** 9),4.2) + 0.299 * EOS.BM3(density_n1[n],4290,309.9 * (10 ** 9),4.59)
                            break
                        elif nth_p_n < 120 * (10 ** 9):
                            pressure_n1[n] = EOS.BM3(density_n1[n], 4101,256 * (10 ** 9),4)
                            break
                        else:
                            pressure_n1[n] = EOS.Vinet(density_n1[n], 4058,221 * (10 ** 9),4.2)
                            break
                        
                    elif list(insert_dict.keys()).index(cutoff) == 2:
                    #This "if block" is for calculations within the Ice.
                        
                        #We calculate the n+1th density and pressure using functions from EoS_Bits 
                        density_n1[n] = EOS.IceDensity(nth_p_n)
                        
                        if nth_p_n < 1 * (10 ** 9):
                            pressure_n1[n] = EOS.Murnaghan(density_n1[n], 930, 9.85 * (10 ** 9), 6.6)
                            break
                        elif nth_p_n < 2.1 * (10 ** 9):
                            pressure_n1[n] = EOS.BM3(density_n1[n], 1271, 14.05 * (10 ** 9), 4)
                            break
                        else:
                            pressure_n1[n] = EOS.BM3(density_n1[n], 1456, 14.9 * (10 ** 9), 5.4)
                            break
    
        if all(abs(np.array(list(density_n1.values())[:1000]) - np.array(density_list[:1000]) < discrepancy)):
            return(list(density_n1.values()), list(pressure_n1.values()))
        #If all densities are self-consistent to 10 kg/m^3, return a list of densities and a list of pressures.
    
        else:
            return(Solver(radii_list, insert_dict, density_list = list(density_n1.values()), discrepancy = discrepancy, calls = calls))
        #If densities are not yet self-consistent, reiterate the process.
        
    return([0], [0])
    #If densities do not converge, then return empty lists.


# Let's try creating a planet using ```planetary_dictionary``` and running our Solver to print out our densities and pressures.
# 
# Let's construct a planet that has a radius of 1 Earth radius. It will be 0.3 parts core and 0.5 parts mantle by depth.

# In[8]:


print('Our dictionary:')
print(dct.planetary_dictionary(1, 0.3, 0.5))
print('\n')

print('Our densities and pressures:')
print(Solver([*range(1, int(r_earth), int(r_earth / 1000))], dct.planetary_dictionary(1, 0.3, 0.5)))


# The dictionary has three keys. Each key represents the boundary at which materials change. In other words, our core has a radius of Key 1 meters. Our mantle begins at the edge of our core and extends to Key 2 meters. Our ice layer begins at the edge of the mantle and extends to Key 3 meters.
# 
# Our list of densities and pressures both start at the innermost radius and end with the outermost radius specified in our ```radii_list```.
# 
# Now, let's plot our density and/or pressure as a function of radius.

# In[9]:


def plotter(radii_list, insert_dict, density_or_pressure):
    """
    Plots either the density or pressure profile of a planet as a function of radius.

    Parameters:
    radii_list (list): A list of radial distances from the center of the planet (in meters).
    insert_dict (dict): A dictionary specifying the planet's internal structure, generated by planetary_dictionary().
    density_or_pressure (str): A string, either 'density' or 'pressure', indicating which quantity to plot.

    Displays:
    A matplotlib plot of the selected quantity versus radius.
    """
    if density_or_pressure.lower() == "density":
        to_plot = Solver(radii_list, insert_dict)[0]
    elif density_or_pressure.lower() == "pressure":
        to_plot = Solver(radii_list, insert_dict)[1]
    plt.plot(radii_list, to_plot)
    plt.show()


# Call the ```plotter``` function with three arguments.
# 
# The first argument is our list of radii. In the following example, our list of radii runs from ```1``` to ```1 * 6.371 * (10 ** 6)``` in increments of ```6371``` (a thousand slices for an Earth-sized planet).
# 
# The second argument specifies the planet's profile. We use a planetary dictionary here.
# 
# The third argument is a string specifying either ```'density'``` or ```'pressure'``` to be plotted.

# In[10]:


plotter([*range(1, int(r_earth), int(r_earth / 1000))], dct.planetary_dictionary(1, 0.3, 0.6), 'density')
plotter([*range(1, int(r_earth), int(r_earth / 1000))], dct.planetary_dictionary(1, 0.3, 0.6), 'pressure')


# 
