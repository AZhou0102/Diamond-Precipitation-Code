#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')


# 
# PREM: Preliminary Reference Earth Model
# =======================================
# 
# The Preliminary reference Earth model (PREM) [Dziewonsky1981]_ is a one-dimensional
# model representing the average Earth properties as a function of planetary radius.  The
# model includes the depth, density, seismic velocities, attenuation (Q) and anisotropic
# parameter ($\eta$) on the boundaries of several Earth layers.  The data is loaded
# into :class:`pandas.DataFrame` objects, which can be used to plot and make computations.
# 
# 

# In[88]:


import rockhound as rh
import matplotlib.pyplot as plt

# Load PREM into a DataFrame
prem = rh.fetch_prem()
print(prem)


# Plot density and velocities
fig, axes = plt.subplots(1, 2, figsize=(9, 5), sharey=True)
fig.suptitle("PREM: Preliminary Reference Earth Model")
ax = axes[0]
prem.plot("density", "depth", legend=False, ax=ax)
ax.invert_yaxis()
ax.set_xlabel("Density [g/cm³]")
ax.set_ylabel("Depth [km]")
ax.grid()
ax = axes[1]
for velocity in ["Vpv", "Vph", "Vsv", "Vsh"]:
    prem.plot(velocity, "depth", legend=False, ax=ax, label=velocity)
ax.grid()
ax.legend()
ax.set_xlabel("Velocity [km/s]")
plt.show()


# In[96]:


get_ipython().run_line_magic('run', 'solve_bm3.ipynb')

import rockhound as rh
import matplotlib.pyplot as plt

# Load PREM into a DataFrame
prem = rh.fetch_prem()



xval_list = prem.loc[:,'depth']
sio2_yval_list = []
mgsio3_yval_list = []
for n in range(len(xval_list)):
    #print(n, prem.loc[:,'density'][n] * 1000)
    add_sio2 = bm3(prem.loc[:,'density'][n] * 1000, 2143.5, 309900000000, 4.59)
    sio2_yval_list.append(add_sio2)

    add_mgsio3 = bm3(prem.loc[:,'density'][n] * 1000, 4101, 256000000000, 4)
    mgsio3_yval_list.append(add_mgsio3)
    
    #add_mg2sio3 = bm3(prem.loc[:,'density'][n] * 1000, 2143.5, 309900000000, 4.59)
    #mg2sio3_yval_list.append(add_mg2sio3)

ax.set_xlabel("Depth [km] (not radius)")
plt.plot(xval_list, sio2_yval_list)
plt.plot(xval_list, mgsio3_yval_list)
plt.show()

