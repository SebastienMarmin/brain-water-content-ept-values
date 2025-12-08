import numpy as np

from ep_water_fitting import EPWaterFitting
from ep_water_uncertainty_MC import WC
from ep_water_fitting_MC import freq_min, freq_max, frequencies, n_lines, norm




def conductivity(age,tissue, frequency):
    wc = WC(age,tissue)
    return EPWaterFitting(frequency, wc)[1]

def permittivity(age,tissue, frequency):
    wc = WC(age,tissue)
    return EPWaterFitting(frequency, wc)[0]


for i,freq in enumerate(frequencies):

    age = 1#np.linspace(0.6,1)
    sigma = conductivity(age,'WM', freq)
    epsilon = permittivity(age,'WM', freq)
    print(sigma)




## Plot
import matplotlib as mpl
import matplotlib.pyplot as plt

print('WAAFF')
print(frequencies)

def myplot(func,name="",unit="-",cmap="binary"):
    fig, axs = plt.subplots(1, 2, figsize=(8, 3),sharey=True,sharex=True)
    Zt = []
    for i, tissue in enumerate(tissues):
        X, Y = np.meshgrid(ages, frequencies)
        print(f"wTF{tissue}")
        print(ages)
        print(tissue)
        print(frequencies[0])
        print(func(ages,tissue,frequencies[0]))

        Z = np.zeros((len(frequencies),len(ages)))

        for l, age in enumerate(ages):
            for k, f in enumerate(frequencies):
                Z[k,l] = func(age,tissue,f)
        Zt += [Z]
    mZ = np.min(tuple(np.min(Z) for Z in Zt))
    MZ = np.max(tuple(np.max(Z) for Z in Zt))
    for i, tissue in enumerate(tissues):
        Z = Zt[i]
        # Plot contour plot with common color scale
        PC = axs[i].pcolor(X, Y/10**6, Z,vmin=mZ,vmax=MZ,cmap=cmap)
        CSF = axs[i].contour(X, Y/10**6, Z,colors="grey",levels= 5)
        
        if name=="Conductivity" and tissue!="CSF":
            fmt = mpl.ticker.StrMethodFormatter("\n{x:1.2f}")
            cl = axs[i].clabel(CSF, inline=False, fmt=fmt, fontsize=6, colors="black")
        else:
            clabels = axs[i].clabel(CSF, fontsize=6, colors="black")
        #im = axs[i].contourf(X, Y, Z, cmap='viridis',vmin=mZ,vmax=MZ)
        axs[i].set_title(tissue)
        axs[i].set_xlabel('Age (y)')
        if i==0:
            axs[i].set_ylabel('Frequency (MHz)')

    

    # Add colorbar
    # bounds = [-1, 2, 5, 7, 12, 15]
    # cmap = mpl.cm.viridis
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(PC, cax=cbar_ax,label='a colorbar label')
    #cbar.add_lines(CSF)
    #cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(mZ, MZ), cmap='viridis'),
    #         ax=axs, orientation='vertical', label='a colorbar label')
    cbar.set_label(name+" ("+unit+")")
    


    
ages = np.linspace(.75,56,100)
tissues = ("WM","GM")#,"CSF")
frequencies = np.linspace(50,350,120)*10**6

myplot(conductivity,"Conductivity", "S/m",cmap="rainbow")
plt.savefig("../article/conductivity.pdf")
plt.close()
myplot(permittivity,"Relative permittivity",cmap="rainbow")
plt.savefig("../article/permittivity.pdf")
    

ages = np.linspace(.75,15,100)
tissues = ("WM","GM")#,"CSF")
frequencies = np.linspace(50,75,120)*10**6

myplot(conductivity,"Conductivity", "S/m",cmap="rainbow")
plt.savefig("../article/conductivity_young.pdf")
plt.close()
myplot(permittivity,"Relative permittivity",cmap="rainbow")
plt.savefig("../article/permittivity_young.pdf")
