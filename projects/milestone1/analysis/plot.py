import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
from astropy.constants import c  
import os 
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import pandas as pd 

import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning )



#   rc and plot params
TICKLABELSIZE = 20
LABELSIZE     = 25
LEGENDSIZE    = 20
TITLESIZE     = 30
LABELPAD      = 1
TITLEPAD      = 15


#   Set rc-params
plt.rc("legend", fontsize=LEGENDSIZE, fancybox=True, loc="best", frameon=True, edgecolor="black")
plt.rc("font", size=25)
plt.rc("axes", titlesize=TITLESIZE, labelsize=LABELSIZE, labelpad=LABELPAD, titlepad=TITLEPAD)
plt.rc("xtick", labelsize=TICKLABELSIZE)
plt.rc("ytick", labelsize=TICKLABELSIZE)
# plt.rc("ytick.minor", pad=1)
plt.rc("lines", linewidth=2.5)
plt.rc('text', usetex=True)
# plt.rc("label", fontsize=20)

plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = 'Times New Roman'



# Paths 
project_path = "/home/vetle/Documents/master_studies/subjects/V23/AST5220/projects/milestone1"


here = os.path.abspath(".")
data_path = project_path + "/data/"
latex_path = project_path + "/m1tex/tables/"
# temp_path = here + "/../../output/plots/temp/"
# fig_path = here +"/figures/"
pdf_path = project_path + "/analysis/figures/"
png_path = pdf_path + "temp/"

def save_push(fig, pdf_name, save=True, push=False, show=False, tight=True, temp=False):
    """
    This function handles whether you want to show,
    save and/or push the file to git.
    Args:
        fig (matplotlib.figure): Figure you want to handle
        pdfname (string): Name of output pdf file
        args (argparse)
    """
    if tight:
        fig.tight_layout()

    if temp:
        fig_name = pdf_name.replace(".pdf", ".png")
        fig_path = png_path
    else:
        fig_name = pdf_name 
        fig_path = pdf_path

    file = fig_path + fig_name

    if save:
        print(f'Saving plot: {fig_name}')
        if not temp:
            fig.savefig(file)
        else:
            fig.savefig(file, dpi=50)
    else:
        plt.show()
    if push:
        os.system(f"git add {file}")
        os.system("git commit -m 'upload plot'")
        os.system("git push")

    plt.close()


def set_ax_info(ax, xlabel, ylabel=False, title=None, legend=True, legendloc='best', legend_size=False):
    """Write title and labels on an axis with the correct fontsizes.
    Args:
        ax (matplotlib.axis): the axis on which to display information
        title (str): the desired title on the axis
        xlabel (str): the desired lab on the x-axis
        ylabel (str): the desired lab on the y-axis
    """
    ax.set_xlabel(xlabel)
    if ylabel != False:
        ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    # ax.tick_params(axis='both', which='major', labelsize=15)
    # ax.yaxis.get_offset_text().set_fontsize(15)
    # try:
        # ax.ticklabel_format(style=style)
    # except AttributeError:
        # pass
    if legend:
        if legend_size:
            ax.legend(loc=legendloc, fontsize=legend_size)
        else:

            ax.legend(loc=legendloc)




# -----------------------------------------------------------------------------
#   General plot code above
#   XXX
#   Specific plot code below
# -----------------------------------------------------------------------------



def load(file, folder=data_path, skiprows=0):
    # Tired of repeating unpack, delimiter, skiprows for all... 
    return np.loadtxt(folder + file, unpack=True, skiprows=skiprows)



def plot_single_param(x, quantity, fname, 
                      mr_eq=None, mL_eq=None, acc=None, 
                      xlabel=None, ylabel=None, 
                      xlim=None, ylim=None,
                      legend=False, legendloc='best', yticks=None,
                      title=None, figsize=(8,6), log=True, 
                      save=True, push=False, temp=False):

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(x, quantity, color='blue')

    if ylim is None:
        ymin = np.min(quantity)
        ymax = np.max(quantity)
    else:
        ymin, ymax = ylim


    if mr_eq is not None:
        ax.vlines(mr_eq, ymin, ymax, ls='dashed', alpha=0.5,
                    color='red', label=r'$\Omega_\mathrm{rel}=\Omega_m$')

    if mL_eq is not None:
        ax.vlines(mL_eq, ymin, ymax, ls='dashed', alpha=0.5,
                    color='green', label=r'$\Omega_m=\Omega_\Lambda$')
    if acc is not None:
        ax.vlines(acc, ymin, ymax, ls='dashed', alpha=0.5, color='purple', label=r'$\ddot{a}=0$')
        
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    if log:
        ax.set_yscale('log')
    
    set_ax_info(ax, xlabel, ylabel, title, legend=legend, legendloc=legendloc)
    if yticks is not None:
        ax.set_yticks(yticks)
    save_push(fig, fname, save, push, temp=temp)
    

def plot_omega_params(x, m, r, L, 
                      fname, title=None,
                      mr_eq=None, ml_eq=None, acc=None,
                      xlabel=r"$x$", ylabel=None,
                      xlim=[-15,3], ylim=[0,1.2], 
                      save=True, push=False, temp=False):

    fig, ax = plt.subplots(figsize=(9,7))

    r_, = ax.plot(x, r, color='red', label=r'$\Omega_\mathrm{r} = \Omega_\gamma + \Omega_\nu$')
    m_, = ax.plot(x, m, color='green', label=r'$\Omega_\mathrm{m}=\Omega_\mathrm{b} + \Omega_\mathrm{CDM}$')
    l_, = ax.plot(x, L, color='orange', label=r'$\Omega_\Lambda$')



    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    leg1 = plt.legend(handles=[r_], loc='upper left' , handlelength=1, fontsize=20)
    leg2 = plt.legend(handles=[m_], loc=(0.44, 0.88)  , handlelength=1, fontsize=20)
    # leg2 = plt.legend(handles=[m_], loc='upper center'  , handlelength=1, fontsize=20)
    leg3 = plt.legend(handles=[l_], loc='upper right', handlelength=1, fontsize=20)

    ax.add_artist(leg1)
    ax.add_artist(leg2)
    ax.add_artist(leg3)
    v_handles = []
    if mr_eq is not None:
        v1 = ax.vlines(mr_eq, *ylim, ls='dashed', alpha=0.5,
                       color='red', label=r'$\Omega_\mathrm{rel}=\Omega_m$')
        v_handles.append(v1)

    if ml_eq is not None:
        v2 = ax.vlines(ml_eq, *ylim, ls='dashed', alpha=0.5,
                       color='green', label=r'$\Omega_m=\Omega_\Lambda$')
        v_handles.append(v2)
        
    if acc is not None:
        v3 = ax.vlines(acc, *ylim, ls='dashed', alpha=0.5, 
                       color='purple', label=r'$\ddot{a}=0$')
        v_handles.append(v3)
        
    if len(v_handles)>0:
        v_leg = plt.legend(handles=v_handles, loc='center left')
        ax.add_artist(v_leg) 

    set_ax_info(ax, xlabel, ylabel, title, legend=False)#, legendloc="center")
    
    save_push(fig, fname, save, push, temp=temp)



def compare_dH_and_ddH_over_H(x, dH_over_H, ddH_over_H, 
                              dH_label, ddH_label, 
                              x_mr_eq, x_mL_eq, x_acc=None,
                              title=None, xlim=[-16,5], ylim=[-1.1, 1.5],
                              save=True, push=False, temp=False):
    
    fig, ax = plt.subplots(figsize=(11,8))
    dH_,  = ax.plot(x, dH_over_H, color='blue', label=dH_label)
    ddH_, = ax.plot(x, ddH_over_H, color='k', label=ddH_label)
    H_handles = [dH_, ddH_]

    l1 = ax.hlines(-1, xmin=xlim[0], xmax=x_mr_eq,  ls='dashed',color='r', label=r'$w=1/3$')
    l2 = ax.hlines(-1/2, xmin=x_mr_eq, xmax=x_mL_eq,ls='dashed',color='green', label=r'$w=0$')
    l3 = ax.hlines(1, xmin=x_mL_eq, xmax=xlim[-1],  ls='dashed',color='orange', label=r'$w=-1$')
    hline_handles = [l1,l2,l3]

    ax.hlines(1, xmin=xlim[0], xmax=x_mr_eq,  ls='dashed',color='r')
    ax.hlines(1/4, xmin=x_mr_eq, xmax=x_mL_eq,ls='dashed',color='green')
    ax.hlines(1, xmin=x_mL_eq, xmax=xlim[-1],  ls='dashed',color='orange')

    l4 = ax.vlines(x_mr_eq, *ylim, ls='dotted', alpha=0.5,
                    color='red', label=r'$\Omega_\mathrm{rel}=\Omega_m$')
    l5 = ax.vlines(x_mL_eq, *ylim, ls='dotted', alpha=0.5,
                    color='green', label=r'$\Omega_m=\Omega_\Lambda$')
    vline_handles = [l4,l5]
    if x_acc is not None:
        l6 = ax.vlines(x_acc, *ylim, ls='dotted', alpha=0.5,
                       color='purple', label=r'$\ddot{a}=0$')
        vline_handles.append(l6)
    
    leg1 = plt.legend(handles=H_handles, loc='center left', fontsize=40, handlelength=1)
    leg2 = plt.legend(handles=hline_handles, loc='lower right')
    leg3 = plt.legend(handles=vline_handles, loc=(0.435,0.8))#, bbox_to_anchor=(-0.1,1))
    ax.add_artist(leg1)
    ax.add_artist(leg2)
    ax.add_artist(leg3)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    set_ax_info(ax, xlabel='$x$', title=title, legend=False)
    # ax.set_xticks([-15,-12,-9, -6, -3, 0, 3])
    ax.set_xticks([-15,-10,-5, 0, 5])


    save_push(fig, 'dH_and_ddH_over_H.pdf', save, push=push, temp=temp)


        
def compare_ddH_over_H(x, ddH_H, H_label, x_mr_eq, x_mL_eq, 
                        title=None, xlim=[-20,5], save=True, push=False, temp=False):
    
    fig, ax = plt.subplots(figsize=(12,10))

    ax.plot(x, ddH_H, label=H_label)
    ax.hlines(1, xmin=xlim[0], xmax=x_mr_eq,  ls='dashed',color='r', label=r'$\omega=1/3$')
    ax.hlines(1/4, xmin=x_mr_eq, xmax=x_mL_eq,ls='dashed',color='green', label=r'$\omega=0$')
    ax.hlines(1, xmin=x_mL_eq, xmax=xlim[-1],  ls='dashed',color='orange', label=r'$\omega=-1$')

    ax.set_xlim(xlim)

    set_ax_info(ax, xlabel='$x$', title=title)
    save_push(fig, 'ddH_over_H.pdf', save, temp=temp)



def plot_t_and_eta(x, t, etac, fname, 
                   mr_eq=None, mL_eq=None, acc=None, 
                   xlim=[-17,5], ylim=[5e-16, 5e2],
                   ylabel=None, xlabel=r"$x$",
                   save=True, temp=False, push=False):

    fig, ax = plt.subplots(figsize=(10,7))

    l1, = ax.plot(x, etac, color='blue', label=r'$\eta(x)/c$')
    l2, = ax.plot(x, t, color='k', label=r'$t(x)$')
    # ylims = ax.get_ylim()
    if mr_eq is not None:
        l3 = ax.vlines(mr_eq, *ylim, ls='dashed', color='red', alpha=0.5, label=r'$\Omega_\mathrm{rad}=\Omega_m$')
    if mL_eq is not None:
        l4 = ax.vlines(mL_eq, *ylim, ls='dashed', color='green', alpha=0.5, label=r'$\Omega_m=\Omega_\Lambda$')
    if acc is not None:
        l5 = ax.vlines(acc, *ylim, ls='dashed', color='purple', alpha=0.5, label=r'$\ddot{a}=0$')

    leg1 = plt.legend(handles=[l1,l2], loc='upper left')
    leg2 = plt.legend(handles=[l3,l4,l5], loc=(0.445, 0.1))
    ax.add_artist(leg1)
    ax.add_artist(leg2)



    # ylabel=r'Time $[\mathrm{Gyr}]$'
    # xlabel=r"$x$"

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_yscale('log')

    set_ax_info(ax, xlabel, ylabel, legend=False)
    save_push(fig, fname, save=save, temp=temp, push=push)



def plot_dL(data, planck=None, fit=[None,None], 
            fname="dL_z_compare.pdf", save=True, temp=False):

    z_data,   dL_data, dL_error = data 
    z_planck, dL_planck         = planck 
    z_fit,    dL_fit            = fit         
    fig, ax = plt.subplots(figsize=(8,7))

    if dL_fit is not None:
        ax.plot(z_fit, dL_fit/z_fit, color='blue', label='MCMC fit')

        ax.plot(z_planck, dL_planck/z_planck, color='green', ls='dashed', alpha=0.7, label='Planck Cosmology')
    else:
        ax.plot(z_planck, dL_planck/z_planck, color='green', label='Planck Cosmology')

    ax.errorbar(z_data, dL_data/z_data, yerr=dL_error/z_data, barsabove=True, fmt='o', 
                capthick=1.5, capsize=5, elinewidth=2, color='r', ms=3,
                label='Data')

    xlabel=r'z'
    ylabel=r'$d_L(\mathrm{z})/\mathrm{z} \:\: [\mathrm{Gpc}]$'

    ax.set_xscale('log')
    ax.set_yscale('log')

    set_ax_info(ax, xlabel, ylabel)#, title="Supernova fit")
    ax.set_ylabel(ylabel, labelpad=10)
    # formatter = ScalarFormatter()
    # formatter.set_scientific(False)
    # ax.yaxis.set_minor_formatter(ScalarFormatter())
    ax.set_yticks([4,5,6,7,8])
    ax.set_yticklabels(['$4$','$5$','$6$','$7$','$8$'])
    # ax.set_axes_locator()
    # ax.yaxis.set_minor_formatter(FormatStrFormatter("%d"))
    # ax.yaxis.set_major_formatter(ScalarFormatter())

    # plt.ticklabel_format(axis='y', style='plain')
    save_push(fig, fname, save=save, temp=temp)



def plot_OmegaM_OmegaLambda_plane(OmegaM, OmegaLambda, chi2_1sigma, chi2_2sigma, chi2_min, 
                                  fname, save=True, temp=False):

    fig, ax = plt.subplots(figsize=(10,8))

    ax.plot(OmegaM[chi2_2sigma], OmegaLambda[chi2_2sigma],'ro',ms=1)
    ax.plot(OmegaM[chi2_1sigma], OmegaLambda[chi2_1sigma],'bo',ms=1)
    ax.plot([], 'ro', label=r'$2\sigma$')
    ax.plot([], 'bo', label=r'$1\sigma$')
    omega_k_zero = np.linspace(0, 1, 20)
    ax.plot(omega_k_zero, 1 - omega_k_zero, 'k--', label=r'$\Omega_k=0$')
    ax.plot(OmegaM[chi2_min], OmegaLambda[chi2_min], 'D', color='orange', ms=10, label='Best fit')

    ax.set_xlim(0,0.8)
    ax.set_ylim(0.1,1.2)
    set_ax_info(ax, xlabel=r"$\Omega_m$", ylabel=r"$\Omega_\Lambda$")
    save_push(fig, fname, save, temp=temp)


def plot_OmegaM_OmegaK_plane(OmegaM, OmegaK, chi2_1sigma, chi2_2sigma, chi2_min, 
                                  fname, save=True, temp=False):

    fig, ax = plt.subplots(figsize=(10,8))

    ax.plot(OmegaM[chi2_2sigma], OmegaK[chi2_2sigma],'ro',ms=1)
    ax.plot(OmegaM[chi2_1sigma], OmegaK[chi2_1sigma],'bo',ms=1)
    ax.plot([], 'ro', label=r'$2\sigma$')
    ax.plot([], 'bo', label=r'$1\sigma$')
    # omega_k_zero = np.linspace(0, 1, 20)
    # ax.plot(omega_k_zero, 1 - omega_k_zero, 'k--', label=r'$\Omega_k=0$')
    ax.plot(OmegaM[chi2_min], OmegaK[chi2_min], 'D', color='orange', ms=10, label='Best fit')

    # ax.set_xlim(0,0.8)
    # ax.set_ylim(0.1,1.2)
    set_ax_info(ax, xlabel=r"$\Omega_m$", ylabel=r"$\Omega_K$")
    save_push(fig, fname, save, temp=temp)

def plot_H0_posterior_pdf(H0, bins, H0_gaussian, fname, bestfit=None, save=True, temp=False):
    fig, ax = plt.subplots(figsize=(12,8))

    ax.plot(bins, H0_gaussian, color='black', label=r"$H_0\sim \mathcal{N}(\mu,\sigma^2)$")
    ax.hist(H0, bins=bins, density=True, color='skyblue', edgecolor='k')
    if bestfit is not None:
        ylim = ax.get_ylim()
        ax.vlines(bestfit, *ylim, ls='solid', color='red', label='Best fit')

    set_ax_info(ax, xlabel=r"$H_0\:[\mathrm{km/s/Mpc}]$")
    save_push(fig, fname, save, temp=temp)

def plot_fit_posterior_pdf(y, bins, fname, mean=None, gaussian=None, save=False, temp=False,
                           param_name=None, xlabel=None):
    
    fig, ax = plt.subplots(figsize=(12,8))
    
    if gaussian is not None:
        ax.plot(bins, gaussian, color='blue', label= param_name + r"$ \sim \mathcal{N}(\mu,\sigma^2)$")
    
    ax.hist(y, bins=bins, density=True, color='green', edgecolor='k')
    
    if mean is not None:
        ylim = ax.get_ylim()
        ax.vlines(mean, *ylim, ls='dashed', color='k', label='Mean')

    set_ax_info(ax, xlabel=xlabel)
    save_push(fig, fname, save, temp=temp)

def bestfit_supernova_table(best_fit, standard_deviations, save=False):
    # Define the variable names, best fit values, and standard deviations
    variables = [r'$\hat{H}_0$', r'$\hat{\Omega}_{m0}$', r'$\hat{\Omega}_{k0}$', r'$\chi_\mathrm{min}^2/N$']
    standard_deviations[0:3] = np.round(standard_deviations[0:3], 4)
    best_fit.append(best_fit.pop(0))
    best_fit = np.round(best_fit, 3)
    
    # Create a dictionary from the variables, best fit values, and standard deviations
    data = {'Parameter': variables, 'Best Fit': best_fit, 'Standard Deviation': standard_deviations}
    # Create a dataframe from the dictionary
    df = pd.DataFrame(data)
    
    # Format the "Best Fit" column to include the standard deviation
    df['Best Fit'] = "$ " + df['Best Fit'].astype(str) + r' \pm ' + df['Standard Deviation'].astype(str) + "$"
    df['Best Fit'][3] = "$" + df['Best Fit'][3].split(" ")[1] + "$"
    df['Best Fit'][0] = df['Best Fit'][0][:-1] + r"\,\mathrm{km/s/Mpc}$" #df['Best Fit'][0].replace()
    
    # Drop the "Standard Deviation" column
    df = df.drop(columns='Standard Deviation')

    # Set table info 
    table_caption   = "The parameter values resulting in the lowest value of $\chi^2$, obtained from fitting to supernova data."
    table_name      = "best_fit_supernova"
    table_label     = "tab:M1:results:" + table_name 
    table_fname     = table_name + ".tex"

    if save:
        buffer = latex_path + table_fname
    else:
        buffer = None 


    # Convert the dataframe to a LaTeX table
    latex_table = df.to_latex(index=False,
                              escape=False,
                              caption=table_caption,
                              label=table_label,
                              position="ht",
                              buf=buffer)
    
    # Print the LaTeX table
    if not save:
        print(latex_table)


def time_table(mr_eq, ml_eq, acc_onset, t0, eta0_over_c, 
               save=False, temp=False, push=False):
    """
    Complete mess...
    """

    # mr_eq = [x_mr_eq, z_mr_eq, t_mr_eq]
    # ml_eq = [x_ml_eq, z_ml_eq, t_ml_eq]
    # acc_onset = [x_acc_onset, z_acc_onset, t_acc_onset]
    # t0 = r"$t_0\:[\mathrm{Gyr}]$" 
    # eta0_row = [eta0_over_c, r"\multicolumn{4}{c}{ }"] 

    row_labels = [r"$x$", 
                  r"$z$", 
                  r"$t\,\mathrm{Gyr}$"]#,

   
    data = np.array([row_labels, mr_eq, acc_onset, ml_eq]).T
    # data = np.vstack([data, t0_row])
   
    col_labels = [r" ", 
                  r"$\Omega_m=\Omega_r$", 
                  r"$\ddot{a}=0$", 
                  r"$\Omega_m=\Omega_\Lambda$"]
    col_fmt = "l|ccc"

    # row_labels = [r"$x$", r"$z$", r"$t\,\mathrm{Gyr}$"]
    # t0_col_name = r"$t_0\:[\mathrm{Gyr}]$"
    # eta_col_name = r"$\eta_0 / c\: [\mathrm{Gyr}]$"


    table_name = "time_valuesNEW"
    table_label = "tab:M1:results:time_values" 
    table_fname = table_name + ".tex"

    table_caption = "Important times during the evolution of the Universe, "
    table_caption += r"expressed in terms of $x$, redshift, $z$ and cosmic time, $t$. "
    table_caption += "In the last two rows we also present today’s time values"

    if save:
        if temp:
            table_fname = "TEMP" + table_fname

        buffer = latex_path + "AVOIDOVERWRITE" + table_fname
    else:
        buffer = None 

    # df = pd.DataFrame(data)
    # df = df.append(pd.Series(['\midrule','\midrule', '\midrule', '\midrule'], index=df.columns), ignore_index=True)
    # print(df)
    # exit()
    table = pd.DataFrame(data).to_latex(
            index=False, 
            header=col_labels, 
            escape=False, 
            column_format=col_fmt,
            multirow=True,
            caption=table_caption,
            label=table_label,
            position="h",
            buf=buffer
        )


    if save:
        print(f"Saving time table for M1:")
        print(f"   {table_fname}")
        
    else:
        print(table)
        print(t0)
        print(eta0_over_c)




    # df = pd.DataFrame(data, index=None)
    # if save:
    #     fname = latex_path+"time_valuesDONTOVERWRITE.tex"
    #     df.style.format("{:.3f}", subset=column_names).hide(axis="index").to_latex(buf=fname, hrules=True)

    # if show:
    #     df.style.format("{:.3f}", subset=column_names).hide(axis="index")
    #     print(df.to_string())



if __name__=='__main__':
    pass 