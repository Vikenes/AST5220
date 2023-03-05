import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
from astropy.constants import c  
import os 
import plot


import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning )

save = False

def data(fname="cosmology.txt"):
    # simulation_data = plot.load("cosmology.txt") # (np.log(1e-10), 5)
    # simulation_data = plot.load("cosmology_times.txt")
    global simulation_data
    global x 
    simulation_data = plot.load(fname)
    x = simulation_data[0]

data()
print(x[0], x[-1])
# x = simulation_data[0]


def load_H_parameters(convert_units=True):

    H_params = simulation_data[1:4] / u.s 

    if convert_units: 
        H_params = H_params.to(100*u.km / u.s / u.Mpc) 

    Hp, dHp_dx, ddHp_ddx = H_params 

    return Hp, dHp_dx, ddHp_ddx


def load_eta_and_t(convert_units=True):
    eta = simulation_data[4] * u.m 
    t   = simulation_data[5] * u.s 

    if convert_units:
        eta = eta.to(u.Mpc)
        t = t.to(u.Gyr)

    return eta, t 


def load_omegas():
    OmegaB, OmegaCDM, OmegaLambda = simulation_data[6:9]
    OmegaR, OmegaNu, OmegaK = simulation_data[9:12]

    OmegaM = OmegaB + OmegaCDM
    OmegaRel = OmegaR + OmegaNu 

    return OmegaM, OmegaRel, OmegaLambda 


def equality_times(idx=False):
    OmegaM, OmegaRel, OmegaLambda = load_omegas()

    Rad_dom_end = np.where(OmegaRel < OmegaM+OmegaLambda)[0][0]
    Lambda_dom_onset = np.where(OmegaLambda > OmegaM+OmegaRel)[0][0]

    MR_eq_idx = np.abs(OmegaRel[:Rad_dom_end+1] - OmegaM[:Rad_dom_end+1]).argmin() 
    ML_eq_idx = np.abs(OmegaM[MR_eq_idx:Lambda_dom_onset+1] - OmegaLambda[MR_eq_idx:Lambda_dom_onset+1]).argmin()
    ML_eq_idx += MR_eq_idx 

        
    if idx:
        return MR_eq_idx, ML_eq_idx
    else:
        return x[MR_eq_idx], x[ML_eq_idx]

def acceleration_onset(idx=False):

    Hp, dHp_dx = load_H_parameters()[0:2]
    a_dot_dot = (np.exp(-x) * Hp * dHp_dx).to((u.Gyr)**(-2)) 

    accleration_onset_idx = np.argmin(np.abs(a_dot_dot))
    if idx:
        return accleration_onset_idx
    else:
        return x[accleration_onset_idx]


def z_of_x(x_):
    return np.exp(-x_) - 1 

def radiation_matter_equality(idx=False):
    OmegaM, OmegaRel, OmegaLambda = load_omegas()
    matter_dom_onset_idx = np.where(OmegaM > OmegaRel+OmegaLambda)[0][0]



def Hp_plot(save):
    ### OK ### 
    Hp = load_H_parameters()[0]
    ylabel = r'$\mathcal{H}\: \left[ \frac{100\,\mathrm{km/s}}{\mathrm{Mpc}} \right] $'

    mr_eq, mL_eq = equality_times()
    print(x[0])
    exit()
    
    plot.plot_single_param(x, Hp, "compare_Hp.pdf", 
                            mr_eq=mr_eq, mL_eq=mL_eq, acc=acceleration_onset(),
                            xlabel=r"$x$", ylabel=ylabel, 
                            xlim=[-16,3], ylim=[1e-1, 1e4], 
                            legend=True, legendloc='lower left', save=save)



def eta_plot(save):
    ### OK ### 
    eta = load_eta_and_t()[0]
    ylabel = r'$\eta\:\:[\mathrm{Mpc}]$'
    mr_eq, mL_eq = equality_times()
    plot.plot_single_param(x, eta, "compare_eta.pdf", 
                            xlabel=r"x", ylabel=ylabel, 
                            xlim=[-15,5], ylim=[0.1, 5e5],
                            mr_eq=mr_eq, mL_eq=mL_eq, log=False,
                            legend=True, save=save)


def eta_H_plot(save):
    ### OK ### 
    Hp = load_H_parameters()[0]
    eta = load_eta_and_t()[0]
    ylabel = r'$\eta \mathcal{H} / c $'

    eta_Hp_over_c = (eta*Hp/c).to(1)
    mr_eq, mL_eq = equality_times()
    plot.plot_single_param(x, eta_Hp_over_c, "compare_eta_H_over_c.pdf", 
                            mr_eq=mr_eq, mL_eq=mL_eq, legend=True,
                            xlabel=r"x", ylabel=ylabel, 
                            xlim=[-16,2], ylim=[0.75, 4], 
                            yticks=[1,2,3,4],
                            log=False, save=save)


def plot_omegas(save):
    ### OK ### 
    OmegaM, OmegaRel, OmegaLambda = load_omegas()
    title = r"$\Omega_i(x)$"

    plot.plot_omega_params(x, OmegaM, OmegaRel, OmegaLambda, xlabel=r'$x$', 
                        fname="omega_i_of_x.pdf", title=title, save=save)
    


def dH_ddH_over_H(save):
    ### OK ### 
    Hp, dHp, ddHp = load_H_parameters()
    m_dom, L_dom = equality_times()

    # dH_label = r"$\frac{1}{\mathcal{H}} \frac{\mathrm{d} \mathcal{H}}{\mathrm{d}x}$"
    # ddH_label = r"$\frac{1}{\mathcal{H}} \frac{\mathrm{d}^2 \mathcal{H}}{\mathrm{d}x^2}$"
    dH_label = r"$\frac{\mathcal{H}'(x)}{\mathcal{H}(x)}$"
    ddH_label = r"$\frac{\mathcal{H}''(x)}{\mathcal{H}(x)}$"


    plot.compare_dH_and_ddH_over_H(x, dHp/Hp, ddHp/Hp, dH_label, ddH_label, m_dom, L_dom, 
                                   title=None, save=save)


def eta_t_plot(save):
    ### OK ### 
    eta, t = load_eta_and_t()
    eta_c = (eta / c).to(u.Gyr)

    mr_eq, mL_eq = equality_times()
    acc = acceleration_onset()

    plot.plot_t_and_eta(x, t, eta_c, fname="t_and_eta_c.pdf", mr_eq=mr_eq, mL_eq=mL_eq, acc=None, save=save)


def luminosity_distance(data):
    
    x_sim =  data[0]
    dL_sim   = (data[-1]*u.m).to(u.Gpc)

    z_sim = z_of_x(x_sim) 
    z_sim_mask = (z_sim <= 1.31) & (z_sim >= 0.005)
    z_sim = z_sim[z_sim_mask]

    dL_sim = dL_sim[z_sim_mask]
    
    return z_sim, dL_sim 


def plot_dL(save, plot_fit=False):
    z, dL, dL_error = plot.load("supernovadata.txt", skiprows=1)
    planck_data     = plot.load("cosmology_dL.txt" , skiprows=1)
    fit_data        = plot.load("bestfit_cosmology_dL.txt" , skiprows=1)

    z_planck, dL_planck = luminosity_distance(planck_data)
    z_fit, dL_fit = luminosity_distance(fit_data)

    data    = [z,        dL, dL_error]
    planck  = [z_planck, dL_planck.value]
    fit     = [z_fit,    dL_fit.value]


    if plot_fit:
        plot.plot_dL(data, planck, fit, fname="dL_z_compare_fitted.pdf", save=save)
    else:
        plot.plot_dL(data, planck, fname='dL_z_compare_planck.pdf', save=save)




def load_supernovafit(burn):
    # supernova_mcmc_results = plot.load("supernovafit_h0_new.txt", skiprows=1)
    supernova_mcmc_results = plot.load("supernovafit.txt", skiprows=1+burn)
    chi2, h, OmegaM, OmegaK = supernova_mcmc_results
    return chi2, h, OmegaM, OmegaK


def supernova_fit_omegas(save, burn=1000):

    chi2, h, OmegaM, OmegaK = load_supernovafit(burn)
   
    OmegaLambda = 1 - OmegaM - OmegaK 
    chi2min = np.min(chi2)
    chi2_1sigma = chi2 < chi2min + 3.53
    chi2_2sigma = chi2 < chi2min + 8.02

    plot.plot_OmegaM_OmegaLambda_plane(OmegaM, OmegaLambda, chi2_1sigma, chi2_2sigma, chi2_min=np.argmin(chi2),
                    fname=f"mcmc_supernova_fit_Nburn{burn}.pdf", save=save)


def supernova_fit_H0_pdf(save, burn=1000):
    chi2, h, OmegaM, OmegaK = load_supernovafit(burn)

    H0 = h * (100*u.km / u.s / u.Mpc)

    H0_mean = np.mean(H0)
    H0_std  = np.std(H0) 
    H0_var  = np.var(H0) 

    print(f"H0 mean: {H0_mean:.5f}")
    print(f"H0 std : {H0_std:.5f}")

    
    H0_min = H0_mean - 4*H0_std 
    H0_max = H0_mean + 4*H0_std 
    N_bins = 100 
    bins   = np.linspace(H0_min, H0_max, N_bins)

    H0_gaussian_distr = (2*np.pi*H0_std**2)**(-1/2) * np.exp(-(bins - H0_mean)**2 / (2 * H0_var))

    plot.plot_H0_posterior_pdf(H0, bins, H0_gaussian_distr,
                                fname=f'H0_pdf_Nburn{burn}.pdf', save=save)


def table():
    data("cosmology_times.txt")
    MR_eq_idx, ML_eq_idx = equality_times(idx=True)
    acc_onset_idx = acceleration_onset(idx=True)
    x_today_idx = np.abs(x).argmin()

    eta, t = load_eta_and_t() 

    Hp, dHp_dx = load_H_parameters()[0:2]



    x_mr = x[MR_eq_idx]
    x_ml = x[ML_eq_idx]
    x_ac = x[acc_onset_idx]

    t_mr = t[MR_eq_idx]
    t_ml = t[ML_eq_idx]
    t_ac = t[acc_onset_idx]

    z_mr = z_of_x(x_mr)
    z_ml = z_of_x(x_ml)
    z_ac = z_of_x(x_ac)

    t_today = t[x_today_idx]
    eta_over_c_today = (eta[x_today_idx] / c).to(u.Gyr) 

    print(f"mr: x={x_mr:.3f}, z={z_mr:.3f}, t={t_mr.to(u.Gyr)/10**(-5):.3f}")
    print(f"ac: x={x_ac:.3f}, z={z_ac:.3f}, t={t_ac:.3f}")
    print(f"ml: x={x_ml:.3f}, z={z_ml:.3f}, t={t_ml:.3f}")
    print(f"t today= {t_today:.3f}")
    print(f"eta0   = {eta_over_c_today:.3f}")

    # print(t[x_today_idx-1:x_today_idx+2])

    mr_eq = [x_mr, z_mr, t_mr.to(u.Gyr).value]
    ml_eq = [x_ml, z_ml, t_ml.value]
    acc   = [x_ac, z_ac, t_ac.value]
    print(x[0], x[-1])

    # print('making table')
    # plot.time_table(mr_eq, ml_eq, acc, t_today, eta_over_c_today, show=True, save=False)

    data()



# plot_dL_best()
# table()

# save=True

# dH_ddH_over_H(save)
# Hp_plot(save)
# eta_plot(save)
# eta_H_plot(save)
# eta_t_plot(save)

# plot_omegas(save)

# plot_dL(save)
# plot_dL(save, True)

# supernova_fit_omegas(save)
# supernova_fit_H0_pdf(save)

