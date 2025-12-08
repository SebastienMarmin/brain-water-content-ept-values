'''
@authors: Alessandro Arduino, Sébastien Marmin
@affiliation: Istituto Nazonale di Ricerca Metrologica, Laboratoire national de métrologie et d'essais
@date: 23/10/2024
'''
import matplotlib as mpl
import numpy as np

import pandas as pd

from matplotlib import pyplot as plt

from .ep_water_fitting_MC import GetReferenceEPsRandom, EpsRFitting, EpsREvaluate, SigmaFitting, SigmaEvaluate, EPWaterFittingRandom

from .ep_water_fitting_MC import EPS0, DISPERSION_PARAMETERS, WATER_CONTENT, DispersionModel, EpsRFitting, EpsREvaluate, SigmaFitting, SigmaEvaluate, freq_min, freq_max, frequencies, n_lines, norm


from .water_content import outputfolder

from os.path import join

from .utils import round_according_to_std, round_with_sense

np.random.seed(0)

# Standard uncertainty for WM (from water_content.py)
u_WM = 0.018839765306671773
u_GM = 0.021528326092851217
stdWC = {'WM':u_WM,'GM':u_GM}

############ From water_content.py

def WC_WM(age):
    return 0.169846315171978 * np.exp(-age/3.3119333386673744) + 0.6866975380846019

def WC_GM(age):
    return 0.09767088107112537 * np.exp(-age/3.7509853460924782) +  0.8291993003134994

def WC_CSF(age):
    return 0.8622907639302754*np.ones_like(age)


def WC(age,tissue):
    if tissue == "CSF":
        return WC_CSF(age)
    elif tissue == "WM":
        return WC_WM(age)
    elif tissue == "GM":
        return WC_GM(age)
    else:
        raise ValueError(f'tissue {tissue} not recognized.')





if __name__=="__main__":

    nmc = 10000

    fig, axes = plt.subplots(1,2, figsize=(10,4))

    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=plt.colormaps["plasma"])
    cmap.set_array([])

    
    std_w = stdWC['GM']

    res = dict()
    water = np.linspace(0.6,1)

    for i,freq in enumerate(frequencies):
        res[freq] = dict()



        eps_r_reps = []
        sigma_reps = []
        for rep in range(nmc):
            w = water + std_w*np.random.normal(0,1)
            eps_r, sigma, water_ref, eps_r_ref, sigma_ref = EPWaterFittingRandom(freq, w)
            eps_r_reps += [eps_r]
            sigma_reps += [sigma]
        m_eps_r = np.mean(eps_r_reps,axis=0)
        m_sigma = np.mean(sigma_reps,axis=0)
        u_eps_r = np.std(eps_r_reps,axis=0)
        u_sigma = np.std(sigma_reps,axis=0)
        res[freq]['m_eps_r'] = m_eps_r
        res[freq]['m_sigma'] = m_sigma
        res[freq]['u_eps_r'] = u_eps_r
        res[freq]['u_sigma'] = u_sigma

        axes[0].plot(water*100, u_eps_r, color=cmap.to_rgba(freq), zorder=i+1)
        axes[0].set_title("Relative permittivity standard uncertainty")
        axes[0].set_ylabel("Relative permittivity standard uncertainty (-)")
        axes[0].set_xlabel("Water content (%)")

        axes[1].plot(water*100, u_sigma, color=cmap.to_rgba(freq), label=f"{int(freq/1e6)} MHz", zorder=i+1)
        axes[1].set_title("Electric conductivity standard uncertainty")
        axes[1].set_ylabel("Electric conductivity standard uncertainty (S/m)")
        axes[1].set_xlabel("Water content (%)")

    axes[0].grid()
    axes[1].grid()

    axes[0].set_axisbelow(True)
    axes[1].set_axisbelow(True)

    axes[1].legend()

    plt.tight_layout()
    plt.savefig(join(outputfolder, "ep_water_uncertainty_MC.pdf"), dpi=300)
    plt.close()

    fig, axes = plt.subplots(1,2, figsize=(10,4))
    for i,freq in enumerate(frequencies):
        
        m_eps_r = res[freq]['m_eps_r']
        m_sigma = res[freq]['m_sigma']
        u_eps_r = res[freq]['u_eps_r']
        u_sigma = res[freq]['u_sigma'] 

        axes[0].plot(water*100, 100*u_eps_r/m_eps_r, color=cmap.to_rgba(freq), zorder=i+1)
        axes[0].set_title("Relative permittivity relative standard uncertainty")
        axes[0].set_ylabel("Relative permittivity relative standard uncertainty (%)")
        axes[0].set_xlabel("Water content (%)")

        axes[1].plot(water*100, 100*u_sigma/m_sigma, color=cmap.to_rgba(freq), label=f"{int(freq/1e6)} MHz", zorder=i+1)
        axes[1].set_title("Electric conductivity relative standard uncertainty")
        axes[1].set_ylabel("Conductivity relative standard uncertainty (%)")
        axes[1].set_xlabel("Water content (%)")

    axes[0].grid()
    axes[1].grid()

    axes[0].set_axisbelow(True)
    axes[1].set_axisbelow(True)


    plt.tight_layout()
    plt.savefig(join(outputfolder, "ep_water_runcertainty_MC.pdf"), dpi=300)
    plt.close()



    def ep_age_uncertainty(age):

        frequencies = (64e6,128e6,230e6,298e6,700e6)
        
        tissues = ('WM','GM')

        max_perm = 0
        max_cond = 0
        min_perm = np.inf
        min_cond = np.inf
        


        fig, axes = plt.subplots(len(frequencies),2, figsize=(10,15),sharex=True)

        cmap = mpl.cm.ScalarMappable(norm=norm, cmap=plt.colormaps["plasma"])
        cmap.set_array([])

        for tissue in tissues:
            for i,freq in enumerate(frequencies):
                lty = 'solid' if tissue == 'WM' else 'dashed'
                mean_w = WC(age,tissue)
                std_w = stdWC[tissue]*np.ones_like(mean_w)
                
                eps_r_reps = []
                sigma_reps = []
                for rep in range(nmc):
                    w = mean_w + std_w*np.random.normal(0,1)
                    assert (np.max(w) <= 1)
                    eps_r, sigma, water_ref, eps_r_ref, sigma_ref = EPWaterFittingRandom(freq, w)
                    eps_r_reps += [eps_r]
                    sigma_reps += [sigma]
                
                eps_r_mean = np.mean(eps_r_reps,0)
                eps_r_lo = np.quantile(eps_r_reps,0.025,axis=0)
                eps_r_hi = np.quantile(eps_r_reps,0.975,axis=0)
                eps_r_U = (eps_r_hi - eps_r_lo)/2

                sigma_mean = np.mean(sigma_reps,0)
                sigma_lo = np.quantile(sigma_reps,0.025,axis=0)
                sigma_hi = np.quantile(sigma_reps,0.975,axis=0)
                sigma_U = (sigma_hi - sigma_lo)/2
                axes[i][0].plot(age, eps_r_mean, color=tuple(fff*.8 for fff in cmap.to_rgba(freq)),linestyle=lty)
                axes[i][0].fill_between(age, eps_r_hi,eps_r_lo, edgecolor=cmap.to_rgba(freq, alpha=0.4),facecolor=cmap.to_rgba(freq, alpha=0.4), linewidth=0.0)
                axes[i][0].set_title(f"{int(freq/1e6)} MHz")
                axes[i][0].set_ylabel("Relative permittivity (-)")
                axes[i][0].set_xlabel("Age (y)" if i == len(frequencies)-1 else '')

                axes[i][1].plot(age, sigma_mean, color=tuple(fff*.8 for fff in cmap.to_rgba(freq)), label=f"{tissue}",linestyle=lty)
                axes[i][1].fill_between(age, sigma_hi,sigma_lo, edgecolor=cmap.to_rgba(freq,alpha=0.4),facecolor=cmap.to_rgba(freq, alpha=0.4), linewidth=0.0)

                axes[i][1].set_ylabel("Conductivity (S/m)")
                axes[i][1].set_xlabel("Age (y)" if i == len(frequencies)-1 else '')

                print(f'tissu {tissue} -- Freq {freq} ')
                dd =   pd.DataFrame(np.concatenate((age.reshape(-1,1), eps_r_mean.reshape(-1,1), eps_r_U.reshape(-1,1),sigma_mean.reshape(-1,1), sigma_U.reshape(-1,1)),axis=1),
                    columns=("age", "eps_r_mean", "eps_r_U","sigma_mean", "sigma_U"))
                print(dd[dd['age'].isin([0,age[5],age[33],age[-1]])])

                axes[i][0].set_axisbelow(True)
                axes[i][1].set_axisbelow(True)
                axes[i][0].margins(x=0)
                axes[i][1].margins(x=0)
                axes[i][0].set_xlim(0,56)
                axes[i][1].set_xlim(0,56)
                if i==0:
                    axes[i][1].legend()

                if min_perm > np.min(eps_r_mean-3.1*eps_r_U/2):
                    min_perm = np.min(eps_r_mean-3.1*eps_r_U/2)
                if min_cond > np.min(sigma_mean-3.1*sigma_U/2):
                    min_cond = np.min(sigma_mean-3.1*sigma_U/2)
                if max_perm < np.max(eps_r_mean+2*eps_r_U/2):
                    max_perm = np.max(eps_r_mean+2*eps_r_U/2)
                if max_cond < np.max(sigma_mean+2*sigma_U/2):
                    max_cond = np.max(sigma_mean+2*sigma_U/2)                    

        for i,freq in enumerate(frequencies):
            axes[i][0].set_ylim(min_perm,max_perm)
            axes[i][1].set_ylim(min_cond,max_cond)
            axes[i][0].grid()
            axes[i][1].grid()


        plt.savefig(join(outputfolder, f"ep_age_uncertainty_{int(np.max(age))}_MC.pdf"), dpi=300)


    ep_age_uncertainty(np.linspace(.75,56,100))



    age=np.array((.75,3,10,56))
    for i,freq in enumerate((64e6,128e6,298e6, 700e6)):
        print(f'freq:{freq}')
        for tissue in ('WM', 'GM'):
            print(f'{tissue}',end='')


            mean_w = WC(age,tissue)
            std_w = stdWC[tissue]*np.ones_like(mean_w)
            
            eps_r_reps = []
            sigma_reps = []
            for rep in range(nmc):
                w = mean_w + std_w*np.random.normal(0,1)
                eps_r, sigma, water_ref, eps_r_ref, sigma_ref = EPWaterFittingRandom(freq, w)
                eps_r_reps += [eps_r]
                sigma_reps += [sigma]
            
            eps_r_mean = np.mean(eps_r_reps,0)
            eps_r_lo = np.quantile(eps_r_reps,0.025,axis=0)
            eps_r_hi = np.quantile(eps_r_reps,0.975,axis=0)
            eps_r_U = (eps_r_hi - eps_r_lo)/2

            sigma_mean = np.mean(sigma_reps,0)
            sigma_lo = np.quantile(sigma_reps,0.025,axis=0)
            sigma_hi = np.quantile(sigma_reps,0.975,axis=0)
            sigma_U = (sigma_hi - sigma_lo)/2

            for i in range(len(eps_r_mean)):
                print('& ', end='')
                print(f'\\makecell[cl]{{$\sigma~$ = {round_according_to_std(sigma_mean[i],sigma_U[i]/2,False)} ± {round_with_sense(sigma_U[i],2)}\\\\ $\\varepsilon_r$ = {round_according_to_std(eps_r_mean[i],eps_r_U[i]/2,False)} ± {round_with_sense(eps_r_U[i],2)} }}')
            print('\\\\')