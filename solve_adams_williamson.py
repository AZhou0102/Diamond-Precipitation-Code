#!/usr/bin/env python
# coding: utf-8

# In[511]:


import numpy as np

r_earth = 6.371 * (10 ** 6)
G = 6.6743015 / (10 ** 11)

def make_density_function(rad, local_densities):
    """
    Turns a list of radii and densities into a dictionary mapping each radius to a density.

    Parameters:
    rad (list): List of radii at which densities are defined.
    local_densities (list): List of densities corresponding to each radius in rad.

    Returns:
    dict: A dictionary mapping each radius to its corresponding density.
    """
    #local_densities should be a list of densities, each corresponding to a single radius.
    density_function = {rad[i]: local_densities[i] for i in range(len(rad))}
    return(density_function)


def find_mass_inside(density_function, current_radius):
    """
    Calculates the mass inside a ball of given radius using a density function.

    Parameters:
    density_function (dict): A dictionary mapping radii to densities.
    current_radius (float): The radius at which to calculate the enclosed mass.

    Returns:
    float: The mass inside the ball of radius current_radius.
    """
    step_size = list(density_function.keys())[1] - list(density_function.keys())[0]
    #Computes the difference between each radius in the radii list.
    mass_inside = 0
    
    for r, local_density in density_function.items():
        if r <= current_radius:
            shell_volume = 4 * np.pi * (r ** 2) * step_size
            shell_mass = shell_volume * local_density
            mass_inside += shell_mass
            #Calculates the mass of each spherical shell within the ball and adds it to the total mass.
            
    return(mass_inside)


def adams_williamson(rad, local_densities):
    """
    Returns a list of gravities and a list of pressures corresponding to each radius within a planet 
    given a list of radii and densities at each radius.

    Parameters:
    rad (list): A list of radii at which to calculate gravity and pressure.
    local_densities (list): A list of densities corresponding to each radius in rad.

    Returns:
    dict: A dictionary mapping each radius to a list of [gravity, pressure].
    """
    density_function = make_density_function(rad, local_densities)
    step_size = list(density_function.keys())[1] - list(density_function.keys())[0]
    gravities_and_pressures = {}
    pressures = []
    total_pressure = 0

    for r, local_density in density_function.items():
        if r <= max(density_function.keys()):
            integrand = (G * find_mass_inside(density_function, r) * local_density) / (r ** 2)
            #Computes the integrand found in the Adams-Williamson equation.
            
            differential_pressure = integrand * step_size
            total_pressure += differential_pressure
            
            #In the dictionary (i.e. gravities_and_pressures), each r corresponds to a list of [gravity, pressure]
            gravities_and_pressures[r] = [integrand / local_density, total_pressure]
    
    for r in gravities_and_pressures.keys():
        gravities_and_pressures[r][1] -= gravities_and_pressures[max(density_function.keys())][1]
        gravities_and_pressures[r][1] = -gravities_and_pressures[r][1]
        
    return(gravities_and_pressures)
    #Returns a dictionary with each radius corresponding to a gravity and a pressure;
    #r: [gravity, pressure]

