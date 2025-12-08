'''
@authors: Alessandro Arduino, Sébastien Marmin
@affiliation: Istituto Nazonale di Ricerca Metrologica, Laboratoire national de métrologie et d'essais
@date: 14/10/2024
'''
import matplotlib as mpl
import numpy as np

from os.path import join

from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as curve_fit

from ep_water_fitting_MC import EPS0, DISPERSION_PARAMETERS, WATER_CONTENT, DispersionModel, EpsRFitting, EpsREvaluate, SigmaFitting, SigmaEvaluate, freq_min, freq_max, frequencies, n_lines, norm

from water_content import outputfolder

def GetReferenceEPs(freq):
    '''
    Provides the electric properties of the reference tissues at the given
    frequency
    '''
    water_ref = []
    eps_r_ref = []
    sigma_ref = []
    for tissue in DISPERSION_PARAMETERS.keys():
        eps_r, sigma = DispersionModel(freq, DISPERSION_PARAMETERS[tissue])
        water_ref.append(WATER_CONTENT[tissue])
        eps_r_ref.append(eps_r)
        sigma_ref.append(sigma)
    water_ref = np.array(water_ref)
    eps_r_ref = np.array(eps_r_ref)
    sigma_ref = np.array(sigma_ref)


    return water_ref, eps_r_ref, sigma_ref



def EPWaterFitting(freq, water_content):
    '''
    Provides the electric properties at the given frequency associated with the
    given water content.
    '''
    water_ref, eps_r_ref, sigma_ref = GetReferenceEPs(freq)
    # fit the relative permittivity
    p = EpsRFitting(water_ref, eps_r_ref)
    eps_r = EpsREvaluate(p, water_content)
    # fit the electric conductivity
    c = SigmaFitting(water_ref, sigma_ref)
    sigma = SigmaEvaluate(c, water_content)
    return eps_r, sigma



if __name__=="__main__":

    fig, axes = plt.subplots(1,2, figsize=(10, 4))

    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=plt.colormaps["plasma"])
    cmap.set_array([])

    for i, freq in enumerate(frequencies):
        water_ref, eps_r_ref, sigma_ref = GetReferenceEPs(freq)

        water = np.linspace(0.6,1)
        eps_r, sigma = EPWaterFitting(freq, water)



        axes[0].plot(water*100, eps_r, color=cmap.to_rgba(freq), zorder=i+1)
        axes[0].scatter(water_ref*100, eps_r_ref, color=cmap.to_rgba(freq), zorder=i+1)
        axes[0].set_title("Relative permittivity")
        axes[0].set_ylabel("Relative permittivity (-)")
        axes[0].set_xlabel("Water content (%)")

        axes[1].plot(water*100, sigma, color=cmap.to_rgba(freq), label=f"{int(freq/1e6)} MHz", zorder=i+1)
        axes[1].scatter(water_ref*100, sigma_ref, color=cmap.to_rgba(freq), zorder=i+1)
        axes[1].set_title("Electric conductivity")
        axes[1].set_ylabel("Electric conductivity (S/m)")
        axes[1].set_xlabel("Water content (%)")

    axes[0].grid()
    axes[1].grid()

    axes[0].set_axisbelow(True)
    axes[1].set_axisbelow(True)

    axes[1].legend()

    plt.tight_layout()
    plt.savefig(join(outputfolder, "ep_water_fitting.pdf"), dpi=300)



