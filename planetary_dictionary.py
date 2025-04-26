#!/usr/bin/env python
# coding: utf-8

# In[1]:


r_earth = 6370 * (10 ** 3)

def planetary_dictionary(earth_rads, iron_part, sio2_part):
    """
    Generates a dictionary representing a planet's composition with three layers: core, mantle, and ice.

    Parameters:
    earth_rads (float): Planet radius in Earth radii.
    iron_part (float): Fraction of the planet's radius made up of the iron core.
    sio2_part (float): Fraction of the planet's radius made up of the mantle.

    Returns:
    dict
        Dictionary with transition radii as keys and initial density guesses as values:
        - Key1: Iron core, density ~7678 kg/m^3.
        - Key2: Silicate mantle, density ~3000 kg/m^3.
        - Key3: Ice layer, density ~930 kg/m^3.
    """
    
    ice_part = 1 - iron_part - sio2_part

    radius = earth_rads * r_earth
    
    iron_radius = iron_part * radius
    sio2_radius = sio2_part * radius
    ice_radius = ice_part * radius

    key1 = int(iron_radius)
    key2 = int(iron_radius) + int(sio2_radius)
    key3 = int(radius)
    
    planetary_dict = {key1: [7678], key2: [3000], key3: [930]}

    #The core resides within the radius of key1.
    #The mantle resides between the radii of key1 and key2.
    #The ice layer resides between the radii of key 2 and key 3.
    #Values for each key represent preliminary density estimates within each layer.
    
    return(planetary_dict)

#Refer to the documentation in Looped_Solver.ipynb for more information.

