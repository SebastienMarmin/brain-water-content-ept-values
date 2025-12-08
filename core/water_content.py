import pandas as pd
import numpy as np
from os.path import join
from pathlib import Path
from matplotlib import colors
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


outputfolder = 'figures'

if __name__ == "__main __":

    np.random.seed(0)

    nmc1 = 10000 # nb of Monte Carlo samples for bootstrap



    Path(outputfolder).mkdir(exist_ok=True)

    file_path = 'alldata.csv'
    data = pd.read_csv(file_path)


    def lighten_color(color, amount=0.5):
        try:
            c = colors.cnames[color]
        except:
            c = color
        c = colors.to_rgb(c)
        c = tuple([(c[i] + (1.0 - c[i]) * amount) for i in range(2)]+[c[-1]])
        return colors.rgb2hex(c)


    def f(x, y0, a, alpha):
        return a * np.exp(-x/alpha) + y0

    parameter_names = ("y0", "a", "alpha")


    def forward(f, x, y, p0, σ):
        popt, pcov = curve_fit(f, x, y, p0, sigma=σ,
                bounds=(np.array([0.0001, 0.0001, 0.0001]), np.inf),
                absolute_sigma=True,maxfev=60000)
        predy = f(x,*popt)
        return y-predy

    def montecarlo(nsample,f,x,y,p0,σ,low,up):
        res = np.empty((nsample,len(x)))
        for i in range(nsample):
            ry = np.random.uniform(y-low,y+up,len(y))
            res[i,:] = forward(f,x,ry,p0,σ)
        return np.std(res)





    for tissue, color in zip(['WM', 'GM'], ['blue', 'red']):
        for group, marker in zip(['control'], ['o']):
            fhere = f
            if tissue == 'GM':
                color = lighten_color(color, 0.3)
            mask = (data["tissu_type"] == tissue) & (data["group"] == group)
            
            low = (data['median_T1'].values[mask]-data['variability1_T1'].values[mask])
            up = (data['variability2_T1'].values[mask]-data['median_T1'].values[mask])
            
            p0 = 0.70, 5, 0.2
            σ = (data['variability2_T1'].values[mask]-data['variability1_T1'].values[mask])/2
            popt, pcov = curve_fit(fhere, data['age'][mask], data['median_T1'][mask], p0, sigma=σ, bounds=(np.array([0.0001, 0.0001, 0.0001]),np.inf),absolute_sigma=True)
            if tissue == 'WM':
                popt_WM = popt
            predx = np.linspace(0.75,56,100)
            predy = fhere(predx, *popt)

            plt.errorbar(data['age'][mask], data['median_T1'][mask],[low, up], marker=marker, color=color, elinewidth=2, capsize=5, alpha=0.7, label=f'{tissue}', linestyle='')
            plt.plot(predx,predy,color=color)
            
            
    plt.xlabel('Age (y)')
    plt.ylabel('T₁ (ms)')
    plt.title('')
    plt.legend()
    plt.grid(True)
    plt.savefig(join(outputfolder, 'T1.pdf'))
    plt.close()

    minage = np.min(data['age'][data['group'] == 'control'])
    maxage = np.max(data['age'][data['group'] == 'control'])


    for tissue, color in zip(['WM', 'GM'], ['blue', 'red']):
        for group, marker in zip(['control'], ['o']):
            fhere = f
            if tissue == 'GM':
                color = lighten_color(color, 0.3)
            mask = (data["tissu_type"] == tissue) & (data["group"] == group)
            print(f'{tissue} - {group}')
            print(f'---------------------')
            x = data['age'][mask]
            y = data['median_WC'][mask]
            σ = (data['variability2_WC'].values[mask]-data['variability1_WC'].values[mask])/2
            p0 = 0.70, 5, 0.2
            popt, pcov = curve_fit(fhere, x, y, p0, sigma=σ, bounds=(np.array([0.0001, 0.0001, 0.0001]),np.inf), absolute_sigma=True)
            
            sigcal = data['sigcal'].values[mask] 
            low = (data['median_WC'].values[mask]-data['variability1_WC'].values[mask])
            up = (data['variability2_WC'].values[mask]-data['median_WC'].values[mask])
            sig_intra = (data['variability2_WC'].values[mask]-data['variability1_WC'].values[mask])/4
            sig_tot = np.sqrt(sigcal**2+sig_intra**2)
            u_tot = montecarlo(nmc1, fhere, x, y, p0, σ, 2*sig_tot, 2*sig_tot)
            print(f'Standard uncertainty for {tissue}-{group}')
            print(u_tot)
            _, caps, _ = plt.errorbar(data['age'][mask], data['median_WC'][mask], [2*sig_tot, 2*sig_tot], marker=marker, color=color, elinewidth=2, capsize=5, alpha=0.7, label=f'{tissue}', linestyle='')
            predx = np.linspace(0.75, 56, 100)
            predy = fhere(predx, *popt)
            plt.fill_between(predx, (predy-2*u_tot), (predy+2*u_tot), alpha=0.2, color=color)
            plt.plot(predx, predy, color=color)
            print("------")
            print(f"WC {tissue} - {group}")
            for pn, v in zip(parameter_names, popt):
                print(f"{pn} = {v}", end="  _  ")
            print(f"fAt {minage} yo = {fhere(minage,*popt)}")
            


    plt.xlabel('Age (y)')
    plt.ylabel('Water content (-)')
    plt.title('')
    plt.legend()
    plt.grid(True)
    plt.savefig(join(outputfolder, 'WC.pdf'))
    plt.close()

    count=0
    for tissue, color in zip(['WM', 'GM'], ['blue', 'red']):
        count+=1
        plt.subplot(2, 1, count)
        for group, marker in zip(['control'], ['o']):
            fhere = f
            if tissue == 'GM':
                color = lighten_color(color, 0.3)
            mask = (data["tissu_type"] == tissue) & (data["group"] == group)
            x = data['age'][mask]
            y = data['median_WC'][mask]
            σ = (data['variability2_WC'].values[mask]-data['variability1_WC'].values[mask])/2
            p0 = 0.70, 5, 0.2
            popt, pcov = curve_fit(fhere, x, y, p0, sigma=σ, bounds=(np.array([0.0001, 0.0001, 0.0001]),np.inf), absolute_sigma=True)
            sigcal = data['sigcal'].values[mask] 
            low = (data['median_WC'].values[mask]-data['variability1_WC'].values[mask])
            up = (data['variability2_WC'].values[mask]-data['median_WC'].values[mask])
            sig_intra = (data['variability2_WC'].values[mask]-data['variability1_WC'].values[mask])/4
            sig_tot = np.sqrt(sigcal**2+sig_intra**2)
            u_tot = montecarlo(nmc1, fhere, x, y, p0, σ, 2*sig_tot, 2*sig_tot)
            preddata = fhere(data['age'][mask],*popt)
            plt.errorbar(data['age'][mask], data['median_WC'][mask]-preddata,[2*sig_tot, 2*sig_tot], marker=marker, color=color, elinewidth=2, capsize=5, alpha=0.7, label=f'{tissue}', linestyle='')
            predx = np.linspace(0.75,56,100)
            predy = fhere(predx,*popt)-fhere(predx,*popt)
            plt.fill_between(predx,(predy-2*u_tot),(predy+2*u_tot), alpha=0.2,color=color)
            plt.plot(predx,predy,color=color)
            plt.ylabel('Residual water content (-)')
            plt.grid(True)
            plt.legend()
    plt.xlabel('Age (y)')
    plt.title('')


    plt.savefig(join(outputfolder, 'WC_residuals.pdf'))
    plt.close()




    print(f'maxage={maxage}')