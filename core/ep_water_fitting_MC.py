'''
@authors: Alessandro Arduino, Sébastien Marmin
@affiliation: Istituto Nazonale di Ricerca Metrologica, Laboratoire national de métrologie et d'essais
@date: 14/10/2024
'''
from os.path import join

import matplotlib as mpl
import numpy as np

from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as curve_fit

from .water_content import outputfolder

import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")



U_LEVEL = 0.1


EPS0 = 8.8541878128e-12
DISPERSION_PARAMETERS = {
    "WM": {
        "eps_inf": 4,
        "delta1": 32,
        "tau1": 7.958e-12,
        "alpha1": 0.1,
        "delta2": 100,
        "tau2": 7.958e-9,
        "alpha2": 0.1,
        "delta3": 40000,
        "tau3": 53.052e-6,
        "alpha3": 0.3,
        "delta4": 35000000,
        "tau4": 7.958e-3,
        "alpha4": 0.02,
        "sigma": 0.02,
    },
    "GM": {
        "eps_inf": 4,
        "delta1": 45,
        "tau1": 7.958e-12,
        "alpha1": 0.1,
        "delta2": 400,
        "tau2": 15.915e-9,
        "alpha2": 0.15,
        "delta3": 200000,
        "tau3": 106.103e-6,
        "alpha3": 0.22,
        "delta4": 45000000,
        "tau4": 5.305e-3,
        "alpha4": 0,
        "sigma": 0.02,
    },
    "CSF": {
        "eps_inf": 4,
        "delta1": 65,
        "tau1": 7.958e-12,
        "alpha1": 0.1,
        "delta2": 40,
        "tau2": 1.592e-9,
        "alpha2": 0,
        "delta3": 0,
        "tau3": 159.155e-6,
        "alpha3": 0,
        "delta4": 0,
        "tau4": 15.915e-3,
        "alpha4": 0,
        "sigma": 2,
    },
}

WATER_CONTENT = {
    "WM": .6949,#0.6574374649832174,#0.68670,
    "GM": .8155,#0.8454231559836316,#.82920,
    "CSF": .9937 #0.9888132510408092
}
# In August 2025, a new database of the IT'IS Fundation has been released. In the new database, besides updated values of the electrical properties that can be used for the calibration, there are also values for the tissue water content. They indicate:
WATER_CONTENT_U = {
    "WM": 2.55/100,
    "GM": 4.36/100,
    "CSF": 2.38/100
}



def DispersionModel(freq, dispersion_parameters):
    '''
    Provides the electric properties at the given frequency associated with the
    given dispersion parameters, according to fourth order Gabriel model.
    '''
    omega = 2*np.pi*freq
    eps_c = dispersion_parameters["eps_inf"]
    for i in range(1, 5):
        eps_c += dispersion_parameters[f"delta{i}"] / \
            (1.0 + (1j*omega*dispersion_parameters[f"tau{i}"])** \
             (1.0 - dispersion_parameters[f"alpha{i}"]))
    eps_c += dispersion_parameters["sigma"] / (1j*omega*EPS0)
    eps_r = np.real(eps_c)
    sigma = -np.imag(eps_c)*omega*EPS0
    return eps_r, sigma

def GetReferenceEPsRandom(freq):
    '''
    Provides the electric properties of the reference tissues at the given
    frequency
    '''
    water_ref = []
    eps_r_ref = []
    sigma_ref = []
    for tissue in DISPERSION_PARAMETERS.keys():
        eps_r, sigma = DispersionModel(freq, DISPERSION_PARAMETERS[tissue])
        w_centre = WATER_CONTENT[tissue]
        w_U = WATER_CONTENT_U[tissue]
        w_min = 0#np.clip(w_centre-w_U,0,1)
        w_max = 1#np.clip(w_centre+w_U,0,1)
        wcandidate = 10
        while wcandidate > w_max or wcandidate < w_min:
            wcandidate = np.random.normal(w_centre,w_U/2)
        water_ref.append(wcandidate)
        eps_r_ref.append(eps_r*(1+np.random.uniform(-U_LEVEL,+U_LEVEL)))
        sigma_ref.append(sigma*(1+np.random.uniform(-U_LEVEL,+U_LEVEL)))
    water_ref = np.array(water_ref)
    eps_r_ref = np.array(eps_r_ref)
    sigma_ref = np.array(sigma_ref)
    return water_ref, eps_r_ref, sigma_ref

def EpsRFitting(water_ref, eps_r_ref):
    '''
    Fits the relative permittivity: e(w) = p1*w**2 + p2*w + p3
    '''
    p = np.polyfit(water_ref, eps_r_ref, 2)
    return p

def EpsREvaluate(p, water_content):
    '''
    Evaluates the relative permittivity: e(w) = p1*w**2 + p2*w + p3
    '''
    eps_r = np.polyval(p, water_content)
    return eps_r

def SigmaFitting(water_ref, sigma_ref):
    '''
    Fits the electric conductivity: c1 + c2*exp(c3*w)
    '''
    f = lambda w, c1,c2,c3: c1 + c2*np.exp(c3*w)
    c, _ = curve_fit(f, water_ref, sigma_ref,maxfev=16000)
    return c

def SigmaEvaluate(c, water_content):
    '''
    Evaluates the electric conductivity: c1*exp(-c2*w) + c3*w + c4
    '''
    f = lambda w, c1,c2,c3: c1 + c2*np.exp(c3*w)
    sigma = f(water_content, *c)
    return sigma


def EPWaterFittingRandom(freq, water_content):
    '''
    Provides the electric properties at the given frequency associated with the
    given water content.
    '''
    nbit = 0
    while nbit < 10:
        try:
            water_ref, eps_r_ref, sigma_ref = GetReferenceEPsRandom(freq)
            # fit the relative permittivity
            p = EpsRFitting(water_ref, eps_r_ref)
            eps_r = EpsREvaluate(p, water_content)
            # fit the electric conductivity
            c = SigmaFitting(water_ref, sigma_ref)
            break
        except RuntimeError:
            nbit += 1
            if nbit == 9: # Virtually never happens
                raise ValueError('Fitting impossible.')
    sigma = SigmaEvaluate(c, water_content)
    return eps_r, sigma, water_ref, eps_r_ref, sigma_ref


freq_min = 50e6
freq_max = 700e6
frequencies = np.array([freq_min, 64e6, 100e6, 128e6, 140e6, 175e6, 210e6, 245e6, 270e6, 298e6, 350e6, freq_max])
n_lines = len(frequencies)
norm = mpl.colors.Normalize(vmin=freq_min, vmax=400e6)


