from os.path import join

import numpy as np


from ep_water_fitting import EPWaterFitting
from ep_water_uncertainty_MC import WC
from ep_water_fitting_MC import freq_min, freq_max, frequencies, n_lines, norm


from water_content import outputfolder

def conductivity(age,tissue, frequency):
    wc = WC(age,tissue)
    return EPWaterFitting(frequency, wc)[1]

def permittivity(age,tissue, frequency):
    wc = WC(age,tissue)
    return EPWaterFitting(frequency, wc)[0]


for i,freq in enumerate(frequencies):

    age = 1
    sigma = conductivity(age,'WM', freq)
    epsilon = permittivity(age,'WM', freq)





## Plot
import matplotlib as mpl
import matplotlib.pyplot as plt



def myplot(func,name="",unit="-",cmap="binary"):
    fig, axs = plt.subplots(1, 2, figsize=(8, 3),sharey=True,sharex=True)
    Zt = []
    for i, tissue in enumerate(tissues):
        X, Y = np.meshgrid(ages, frequencies)

        Z = np.zeros((len(frequencies),len(ages)))

        for l, age in enumerate(ages):
            for k, f in enumerate(frequencies):
                Z[k,l] = func(age,tissue,f)
        Zt += [Z]
    mZ = np.min(tuple(np.min(Z) for Z in Zt))
    MZ = np.max(tuple(np.max(Z) for Z in Zt))
    for i, tissue in enumerate(tissues):
        Z = Zt[i]
        PC = axs[i].pcolor(X, Y/10**6, Z,vmin=mZ,vmax=MZ,cmap=cmap)
        cntr = axs[i].contour(X, Y/10**6, Z,colors="grey",levels= 5)
        
        if name=="Conductivity":
            fmt = mpl.ticker.StrMethodFormatter("\n{x:1.2f}")
            cl = axs[i].clabel(cntr, inline=False, fmt=fmt, fontsize=6, colors="black")
        else:
            clabels = axs[i].clabel(cntr, fontsize=6, colors="black")
        axs[i].set_title(tissue)
        axs[i].set_xlabel('Age (y)')
        if i==0:
            axs[i].set_ylabel('Frequency (MHz)')

 
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(PC, cax=cbar_ax,label='a colorbar label')
    cbar.set_label(name+" ("+unit+")")
    


    
ages = np.linspace(.75,56,100)
tissues = ("WM","GM")
frequencies = np.linspace(50,350,120)*10**6

myplot(conductivity,"Conductivity", "S/m",cmap="rainbow")
plt.savefig(join(outputfolder,"conductivity.pdf"))
plt.close()
myplot(permittivity,"Relative permittivity",cmap="rainbow")
plt.savefig(join(outputfolder,"permittivity.pdf"))
    

ages = np.linspace(.75,15,100)
tissues = ("WM","GM")
frequencies = np.linspace(50,75,120)*10**6

myplot(conductivity,"Conductivity", "S/m",cmap="rainbow")
plt.savefig(join(outputfolder,"conductivity_young.pdf"))
plt.close()
myplot(permittivity,"Relative permittivity",cmap="rainbow")
plt.savefig(join(outputfolder,"permittivity_young.pdf"))
