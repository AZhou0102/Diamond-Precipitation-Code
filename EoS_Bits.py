#!/usr/bin/env python
# coding: utf-8

# In[36]:


import numpy as np
from matplotlib import pyplot as plt


# In[37]:


def BM3(rho, rho0, B0, B1):
    """Calculate pressure using the 3rd-order Birch-Murnaghan equation of state.
    
    Parameters:
    rho (float): Density at pressure P (kg/m^3).
    rho0 (float): Initial (zero-pressure) density (kg/m^3).
    B0 (float): Bulk modulus at zero pressure (Pa).
    B1 (float): Pressure derivative of the bulk modulus.

    Returns:
    Pressure (float) corresponding to the given density (Pa)."""
    
    foo = rho/rho0
    return 1.5*B0*(foo**(7/3)-foo**(5/3))*(1+0.75*(B1-4)*(foo**(2/3)-1))

def Vinet(rho, rho0, B0, B1):
    """Calculate pressure using the Vinet equation of state.
    
    Parameters:
    rho (float): Density at pressure P (kg/m^3).
    rho0 (float): Initial (zero-pressure) density (kg/m^3).
    B0 (float): Bulk modulus at zero pressure (Pa).
    B1 (float): Pressure derivative of the bulk modulus.

    Returns:
    Pressure (float) corresponding to the given density (Pa)."""
    
    eta = (rho0/rho)**(1/3)
    return 3*B0*((1-eta)/eta**2)*np.exp(1.5*(B1-1)*(1-eta))

def Murnaghan(rho, rho0, B0, B1):
    """Calculate pressure using the Murnaghan equation of state.
    
    Parameters:
    rho (float): Density at pressure P (kg/m^3).
    rho0 (float): Initial (zero-pressure) density (kg/m^3).
    B0 (float): Bulk modulus at zero pressure (Pa).
    B1 (float): Pressure derivative of the bulk modulus.

    Returns:
    Pressure (float) corresponding to the given density (Pa)."""
    return (B0/B1)*((rho0/rho)**(-1*B1)-1)

def DensityFromP(P, rho0, B0, B1, form, thresh=0.01):
    """Calculate density from pressure using a specified equation of state.

    Parameters:
    P (float): Pressure (Pa).
    rho0 (float): Initial (zero-pressure) density (kg/m^3).
    B0 (float): Zero pressure bulk modulus (Pa).
    B1 (float): Pressure derivative of the bulk modulus.
    form (str): Name of the equation of state ('bm3', 'vinet', or 'murnaghan').
    thresh (float, optional): Convergence threshold for pressure (Pa). Default is 0.01 Pa.

    Returns:
    Estimated density (float) at pressure P (kg/m^3)."""
    # check form is supported and assign EoS variable to call that function
    supported_EoSs = {'bm3':BM3, 'vinet':Vinet, 'murnaghan':Murnaghan}
    assert form.lower() in supported_EoSs.keys()
    EoS = supported_EoSs[form.lower()]
    
    # find rough bounds for density
    rho_upper = 2*rho0  #initialize for loop
    rho_lower = 0.9*rho0  # below rho0 to avoid weirdness at 0 pressure
    while True:
        P_upper = EoS(rho_upper, rho0, B0, B1)
        P_lower = EoS(rho_lower, rho0, B0, B1)
        if P<P_upper and P>P_lower:  # bounds are good
            break
        if P>P_upper:  #highest density too low
            rho_lower = rho_upper
            rho_upper += rho0
        if P<P_lower:  # lowest density to too high
            P_upper = P_lower
            P_lower /= 10
    
    # iterate until P_guess within thresh of P when calculated using rho_guess
    rho_guess = (rho_upper+rho_lower)/2
    while True:
        P_guess = EoS(rho_guess, rho0, B0, B1)
        if np.abs(P-P_guess) < thresh:  # rho_guess is good
            break
        if P_guess > P: # rho_guess is too high
            rho_upper = rho_guess
        if P_guess < P:  # rho_guess is too low
            rho_lower = rho_guess
        
        #  recalculate rho_guess with new bounds
        rho_guess = (rho_upper+rho_lower)/2
    
    return rho_guess    


# In[38]:


# EoS for individual compoonents
def IceDensity(P):
    """Returns the density of ice at a given pressure via corresponding EoS based on phase."""
    
    ice_params = {
        'Ih':{'rho0':930, 'B0':9.85 * (10 ** 9), 'B1':6.6, 'form':'Murnaghan'},
        'VI':{'rho0':1271, 'B0':14.05 * (10 ** 9), 'B1':4, 'form':'BM3'},
        'VII':{'rho0':1456, 'B0':14.9 * (10 ** 9), 'B1':5.4, 'form':'BM3'}
    }
    
    # pick phase according to pressure
    if P<1 * (10 ** 9):
        return DensityFromP(P, **ice_params['Ih'])
    elif P<2.1 * (10 ** 9):
        return DensityFromP(P, **ice_params['VI'])
    else:
        return DensityFromP(P, **ice_params['VII'])
        

def RockDensity(P):
    """Return the density of rock (MgSiO3) at a given pressure.
    Considers phase transitions accounting for a mixture of Mg2SiO4 and SiO2.

    Parameters:
    P (float): Pressure (Pa).

    Returns:
    Density (float) of the rock at the given pressure (kg/m^3)."""
    
    ##############
    # NOTES:
    # Mg:Si is 1:1, meaning that below bridgmanite fomration at 23 GPa will have mix
    # of Mg2SiO4 and SiO2, each with several phase trasnitions.  Mix is equimolar.
    # Molar mass of Mg2SiO4: 140.69 g
    # Molar mass of SiO2: 60.083 g
    # So mass (and density) fraction is 0.299 SiO2 : 0.701
    part_SiO2 = 0.299
    part_Mg2SiO4 = 0.701
    
    rock_params = {
        'forsterite':{'rho0':3221,'B0':125 * (10 ** 9),'B1':4,'form':'BM3'},
        'wadsleyite':{'rho0':3491,'B0':160 * (10 ** 9),'B1':4,'form':'BM3'},
        'ringwoodite':{'rho0':3548,'B0':182 * (10 ** 9),'B1':4.2,'form':'BM3'},
        'quartz':{'rho0':2648,'B0':37.4 * (10 ** 9),'B1':6.2,'form':'BM3'},
        'coesite':{'rho0':2921,'B0':96 * (10 ** 9),'B1':8.4,'form':'BM3'},
        'stichovite':{'rho0':4290,'B0':309.9 * (10 ** 9),'B1':4.59,'form':'BM3'},
        'bridgmanite':{'rho0':4101,'B0':256 * (10 ** 9),'B1':4,'form':'BM3'},
        'ppv':{'rho0':4058,'B0':221 * (10 ** 9),'B1':4.2,'form':'vinet'}
    }
    
    # Mg2SiO4, SiO2, and MgSiO3 all have phase traisitions
    # SiO2: quartz -2.5GPa-> coesite -8GPa-> stichovite
    # Mg2SiO4: forsterite -14GPa-> wadsleyite -18GPa-> ringwoodite
    # MgSiO3 forms at 23 GPa as bridgmanite -120GPa-> PPV
    
    # Pick phases according to pressure
    if P<2.5 * (10 ** 9):  # quartz plus forsterite
        rho_quartz = DensityFromP(P, **rock_params['quartz'])
        rho_forsterite = DensityFromP(P, **rock_params['forsterite'])
        return part_SiO2*rho_quartz + part_Mg2SiO4*rho_forsterite
    elif P<8 * (10 ** 9):  # coesite plus forsterite
        rho_coesite = DensityFromP(P, **rock_params['coesite'])
        rho_forsterite = DensityFromP(P, **rock_params['forsterite'])
        return part_SiO2*rho_coesite + part_Mg2SiO4*rho_forsterite
    elif P<14 * (10 ** 9):  # stichovite plus forsterite
        rho_stichovite = DensityFromP(P, **rock_params['stichovite'])
        rho_forsterite = DensityFromP(P, **rock_params['forsterite'])
        return part_SiO2*rho_stichovite + part_Mg2SiO4*rho_forsterite
    elif P<18 * (10 ** 9):  # stichovite plus wadsleyite
        rho_stichovite = DensityFromP(P, **rock_params['stichovite'])
        rho_wadsleyite = DensityFromP(P, **rock_params['wadsleyite'])
        return part_SiO2*rho_stichovite + part_Mg2SiO4*rho_wadsleyite
    elif P<23 * (10 ** 9):  # stichovite plus ringwoodite
        rho_stichovite = DensityFromP(P, **rock_params['stichovite'])
        rho_ringwoodite = DensityFromP(P, **rock_params['ringwoodite'])
        return part_SiO2*rho_stichovite + part_Mg2SiO4*rho_ringwoodite
    elif P<120 * (10 ** 9):  # bridgmanite
        return DensityFromP(P, **rock_params['bridgmanite'])
    else:   # P>120 and have ppv phase
        return DensityFromP(P, **rock_params['ppv'])
    

def CoreDensity(P):
    """Return the density of Fe93Si7 core alloy at a given pressure using the Vinet EoS (based on Wicks 2018)
    
    Parameters:
    P (float): Pressure (Pa).

    Returns:
    Density (float) of the core alloy at the given pressure (kg/m^3)
    """
    
    core_params = {
        'Fe93Si7':{'rho0':7678, 'B0':136.2 * (10 ** 9), 'B1':5.97, 'form':'vinet'}
    }
    
    return DensityFromP(P, **core_params['Fe93Si7'])



# In[39]:


# make vectorized versions of each
vIceDensity = np.vectorize(IceDensity)
vRockDensity = np.vectorize(RockDensity)
vCoreDesnity = np.vectorize(CoreDensity)


# pressures = np.linspace(0, 300 * (10 ** 9), 600)
# plt.plot(pressures, vIceDensity(pressures), label='ice')
# plt.plot(pressures, vRockDensity(pressures), label='rock')
# plt.plot(pressures, vCoreDesnity(pressures), label='core')
# plt.legend()

# In[ ]:




