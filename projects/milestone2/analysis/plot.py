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
# project_path = os.path.abspath(".")
project_path = "/home/vetle/Documents/master_studies/subjects/V23/AST5220/projects/milestone2"

data_path = project_path + "/data/"
latex_path = project_path + "/report/tables/"
fig_path = project_path + "/analysis/figures/"


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
        fig_name = "TEMP" + pdf_name.replace(".pdf", ".png")
    else: 
        fig_name = pdf_name 
    
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



def plot_quantity_with_derivatives(x, y, dy, ddy,
                                   y_legend, dy_legend, ddy_legend, 
                                   fname, xlabel=r"$x$", ylabel=None, 
                                   xlim=None, ylim=None, 
                                   legendloc='best', yticks=None, 
                                   figsize=(8,6), log=True, 
                                   save=True, temp=False):

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(x, y  ,ls='solid' , color='blue'  , label=y_legend)
    ax.plot(x, dy ,ls='dashed', color='orange', label=dy_legend)
    ax.plot(x, ddy,ls='dotted', color='green' , label=ddy_legend)


    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    if log:
        ax.set_yscale('log')
    
    set_ax_info(ax, xlabel, ylabel, legendloc=legendloc)

    if yticks is not None:
        ax.set_yticks(yticks)

    save_push(fig, fname, save, temp=temp)
    


def time_table(mr_eq, ml_eq, acc_onset, t0, eta0, save=False, show=False):
    """
    Complete mess...
    """



    eta_col_name = r"$\eta_0 / c\: [\mathrm{Gyr}]$"
    t0_col_name = r"$t_0\:[\mathrm{Gyr}]$"
    column_names = [r"$\Omega_m=\Omega_r$", r"$\Omega_m=\Omega_\Lambda$", r"Acceleration onset", t0_col_name, eta_col_name]

    row_labels = np.asarray([r"$x$", r"$z$", r"$t$"])
    data = {
        r"var"                      : row_labels, 
        r"$\Omega_m=\Omega_r$"      : mr_eq,
        r"$\Omega_m=\Omega_\Lambda$": ml_eq,
        r"Acceleration onset"       : acc_onset,
        t0_col_name                 : t0,
        eta_col_name                : eta0
    }


    df = pd.DataFrame(data, index=None)
    if save:
        fname = latex_path+"time_valuesDONTOVERWRITE.tex"
        df.style.format("{:.3f}", subset=column_names).hide(axis="index").to_latex(buf=fname, hrules=True)

    if show:
        df.style.format("{:.3f}", subset=column_names).hide(axis="index")
        print(df.to_string())



if __name__=='__main__':
    pass 