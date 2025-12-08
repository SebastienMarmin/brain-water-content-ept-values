'''
@authors: Alessandro Arduino, Sébastien Marmin
@affiliation: Istituto Nazonale di Ricerca Metrologica, Laboratoire national de métrologie et d'essais
@date: 23/10/2024
'''
import matplotlib as mpl
import numpy as np

import pandas as pd

from matplotlib import pyplot as plt

from core.ep_water_fitting_MC import GetReferenceEPsRandom, EpsRFitting, EpsREvaluate, SigmaFitting, SigmaEvaluate, EPWaterFittingRandom

from core.ep_water_fitting_MC import EPS0, DISPERSION_PARAMETERS, WATER_CONTENT, DispersionModel, EpsRFitting, EpsREvaluate, SigmaFitting, SigmaEvaluate, freq_min, freq_max, frequencies, n_lines, norm


from core.water_content import outputfolder

from os.path import join

from core.utils import round_according_to_std, round_with_sense

from core.ep_water_uncertainty_MC import u_GM, u_WM, stdWC, WC_WM, WC_GM, WC_CSF, WC





if __name__=="__main__":

    print("Enter the frequency in hertz (e.g. ‘128e6’):")
    freq = float(input())
    print("Enter the age in year (e.g. ‘3.5’):")
    age = float(input())

    nmc = 10000

    for tissue in ('WM', 'GM'):
        print('White matter' if tissue == "WM" else 'Grey matter')


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


        print(f'σ  = {round_according_to_std(sigma_mean,sigma_U/2,False)} ± {round_with_sense(sigma_U,2)} S/m')
        print(f'εᵣ = {round_according_to_std(eps_r_mean,eps_r_U/2,False)} ± {round_with_sense(eps_r_U,2)}.')
        print('')