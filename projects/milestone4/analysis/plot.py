import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
from astropy.constants import c  
import os 
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import pandas as pd 

import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning )
warnings.filterwarnings("ignore", category=FutureWarning )

#   rc and plot params
TICKLABELSIZE = 20
LABELSIZE     = 25
LEGENDSIZE    = 20
TITLESIZE     = 30
LABELPAD      = 1
TITLEPAD      = 15


#   Set rc-params
plt.rc("legend", fontsize=LEGENDSIZE, fancybox=True, loc="best", 
       frameon=True, edgecolor="black", framealpha=1)
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
project_path = "/home/vetle/Documents/master_studies/subjects/V23/AST5220/projects/milestone4"

data_path = project_path + "/data/"
latex_path = project_path + "/report/tables/"
pdf_path = project_path + "/analysis/figures/"
png_path = pdf_path + "temp/" 


# Global parameters used in all functions 



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
    pdf_name += ".pdf"

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


def set_ax_info(ax, xlabel, ylabel=False, title=None, legend=True, 
                double_legends=[], legendloc='best', legend_size=False,
                xlim=None, ylim=None, yticks=None,
                ypad=None):
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
        # if len(legends)==2:
            # l1 = legends[0]
            # l2 = legends[1]

        if legend_size:
            ax.legend(loc=legendloc, fontsize=legend_size)
        else:
            ax.legend(loc=legendloc)



    if ylim is not None:
        ax.set_ylim(ylim)
    
    if xlim is not None:
        ax.set_xlim(xlim)

    if ypad is not None:
        ax.set_ylabel(ylabel,labelpad=ypad)

    if yticks is not None:
        ax.set_yticks(yticks)

# -----------------------------------------------------------------------------
#   General plot code above
#   XXX
#   Specific plot code below
# -----------------------------------------------------------------------------



def plot_C_ell(x, y,
               fname, 
               logx=True, 
               logy=True,
               xlabel=r"$\ell$", 
               keq=None,
               ylabel=None, 
               xlim=None, 
               ylim=None, 
               legendloc='best', 
               yticks=None, 
               ypad=None,
               figsize=(10,6), 
               save=True,push=False, temp=False):


    fig, ax = plt.subplots(figsize=figsize)


    
    ax.plot(x, y, color='blue'  )




    
    set_ax_info(ax, xlabel, ylabel, legend=False, 
                ylim=ylim, xlim=xlim, yticks=yticks,
                ypad=ypad)
    if logy:
        ax.set_yscale('log')

    if logx:
        ax.set_xscale('log')



    save_push(fig, fname, save, push, temp=temp)
    



def plot_matter_PS(x, y, k_eq,
                   fname,  
                   logx=True, logy=True, 
                   xlabel=None, ylabel=None,  
                   xlim=None, ylim=None,  
                   ypad=None, 
                   figsize=(10,6),  
                   save=True,push=False, temp=False):


    fig, ax = plt.subplots(figsize=figsize)


    
    ax.plot(x, y, color='blue')

    if ylim is None:
        ylim_vline = ax.get_ylim()
    else:
        ylim_vline = ax.get_ylim()

    if k_eq is not None:
        ax.vlines(k_eq.value, *ylim_vline, colors='red', ls='dashed', label=r"$k_\mathrm{eq}$")



    
    set_ax_info(ax, xlabel, ylabel, legend=True, 
                ylim=ylim, xlim=xlim, ypad=ypad)

    if logy:
        ax.set_yscale('log')

    if logx:
        ax.set_xscale('log')



    save_push(fig, fname, save, push, temp=temp)


def plot_Thetas(x, y,
               ell_values, 
               fname,
               logx=True, 
               logy=True,
               xlabel=r"$k\eta_0$", 
               ylabel=None, 
               xlim=None, 
               ylim=None, 
               legendloc='best', 
               yticks=None, 
               ypad=None,
               figsize=(10,6), 
               save=False,push=False, temp=False):


    fig, ax = plt.subplots(figsize=figsize)
    print(fname)

    for idx, l in enumerate(ell_values):
        ax.plot(x, y[idx], label=fr"$\ell=${str(l)}")




    
    set_ax_info(ax, xlabel, ylabel, legend=True, 
                ylim=ylim, xlim=xlim, yticks=yticks,
                ypad=ypad)
    if logy:
        ax.set_yscale('log')

    if logx:
        ax.set_xscale('log')



    save_push(fig, fname, save, push, temp=temp)


def comp_C_ell(x, y,
               xlabel=r"$\ell$", 
               keq=None,
               ylabel=None, 
               xlim=None, 
               ylim=None, 
               legendloc='best', 
               yticks=None, 
               ypad=None,
               figsize=(10,6), 
               save=True,push=False, temp=False):


    fig, ax = plt.subplots(figsize=figsize)


    
    ax.plot(x, y, color='blue'  )




    
    set_ax_info(ax, xlabel, ylabel, legend=False, 
                ylim=ylim, xlim=xlim, yticks=yticks,
                ypad=ypad)



    # save_push(fig, pdf_name=None, save=save, push=push, temp=temp)
    plt.show()

if __name__=='__main__':
    pass 