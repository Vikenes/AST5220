import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
from astropy.constants import c  
import os 
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import pandas as pd 
from cycler import cycler

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
plt.rc("lines", linewidth=2.0)
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
                xlim=None, ylim=None, yticks=None, xticks=None,
                ypad=None, xpad=None):
    """Write title and labels on an axis with the correct fontsizes.
    Args:
        ax (matplotlib.axis): the axis on which to display information
        title (str): the desired title on the axis
        xlabel (str): the desired lab on the x-axis
        ylabel (str): the desired lab on the y-axis
    """
    ax.set_xlabel(xlabel, labelpad=xpad)
    if ylabel != False:
        ax.set_ylabel(ylabel, labelpad=ypad)
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

    if yticks is not None:
        ax.set_yticks(yticks)

    if xticks is not None:
        ax.set_xticks(xticks)

# -----------------------------------------------------------------------------
#   General plot code above
#   XXX
#   Specific plot code below
# -----------------------------------------------------------------------------



def plot_C_ell(ell, C_ell,
               fname, 
               planck=None,
               peaks=None, troughs=None,
               logx=True, 
               logy=True,
               xlabel=r"$\ell$", 
               keq=None,
               ylabel=None, 
               xlim=None, 
               ylim=None, 
               legendloc='best', 
               yticks=None, 
               xticks=None,
               ypad=None,
               fill=False,
               figsize=(10,6), 
               save=True,push=False, temp=False):


    fig, ax = plt.subplots(figsize=figsize)

    
    if fill:
        T0 = 2.7255
        normfactor = ell*(ell+1)/(2*np.pi) #* (1e6*T0)**2 
        C_ell_true = C_ell #/ normfactor 
        y1 = C_ell - C_ell *(2/(ell*(ell+1)))**(1/2) 
        y2 = C_ell + C_ell *(2/(ell*(ell+1)))**(1/2) 
        ax.fill_between(ell, y1, y2, color='palegreen', alpha=1)
        # plt.xscale('log')
        # plt.show()
        
    if peaks is None:
        cl_color = "blue"
    else:
        cl_color = "black"
    cl, = ax.plot(ell, C_ell, color=cl_color, label='Prediction')
    
    if planck is not None:
        ell_planck, C_ell_planck, error_planck = planck
        err = ax.errorbar(ell_planck, C_ell_planck, error_planck, barsabove=True, fmt='x',
                    capthick=1.5, capsize=3, elinewidth=1, color='red', ms=5, 
                    label='Planck 2018')
        leg1 = plt.legend(handles=[cl, err], loc='upper left')
        ax.add_artist(leg1)


    ylim = ax.get_ylim()
    
    set_ax_info(ax, xlabel, ylabel, legend=False, 
                ylim=ylim, xlim=xlim, yticks=yticks,
                ypad=ypad, legendloc='upper left')
                
    if peaks is not None:
        # n = len(peaks) + len(troughs)
        # new_colors = [plt.get_cmap('jet')(1. * i/n) for i in range(n)] 
        # custom_cycler = (cycler(color=palettes.tue_plot))
        # ax.set_prop_cycle(custom_cycler)

        peak_handles = []
        
        ls_ = ['--', '-.', ':']
        # labels = [rf"$\ell = {{{peak:.0f}}} $" for peak in peaks]
        for idx, peak in enumerate(peaks):
            label = rf"$\ell = {{{peak:.0f}}} $"
            line = ax.vlines(peak, *ylim, ls=ls_[idx], colors='red', label=label)
            peak_handles.append(line)

        # print(lines.get_label())
        # exit()
        # handles = [l for l in line]

        leg_peaks = plt.legend(handles=peak_handles, loc='upper left')
        ax.add_artist(leg_peaks)

    if troughs is not None:

        trough_handles = []
        for idx, trough in enumerate(troughs):
            label = rf"$\ell = {{{trough:.0f}}} $"
            line = ax.vlines(trough, *ylim, ls=ls_[idx], colors='blue', label=label)
            trough_handles.append(line)

        leg_troughs = plt.legend(handles=trough_handles, loc='center left')
        ax.add_artist(leg_troughs)


    if logy:
        ax.set_yscale('log')

    if logx:
        ax.set_xscale('log')
    
    if xticks is not None:
        ticklabels=[f"${tick}$" for tick in xticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(ticklabels)


    
    save_push(fig, fname, save, push, temp=temp)
    

def plot_C_ell_components(ell, C_ell, C_ell_components,
                          ell_planck, C_ell_planck,
                          error_planck, 
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
                          xticks=None,
                          ypad=None,
                          fill=False,
                          figsize=(10,6), 
                          save=True,push=False, temp=False):


    fig, ax = plt.subplots(figsize=figsize)


        

    SW, ISW, Doppler, Quadrupole = C_ell_components

    # C_planck = ax.errorbar(ell_planck, C_ell_planck, error_planck, barsabove=True, fmt='o',
    #             capthick=1.5, capsize=5, elinewidth=2, color='gold', ms=3, alpha=1,
    #             label='Planck 2018')
    scale_1 = 1 
    scale_2 = 1 
    scale_3 = 1 
    scale_4 = 1 
    C4, = ax.plot(ell, Quadrupole * scale_4 , ls='solid', alpha=1, lw=1, color='purple', label="Quadrupole")
    Ctot, = ax.plot(ell, C_ell, ls='dashed', color='blue', alpha=0.7, label='Prediction')
    C1, = ax.plot(ell, SW * scale_1         , ls='solid', color='green', label="SW")
    C2, = ax.plot(ell, ISW * scale_2        , ls='solid', color='red' , label="ISW")
    C3, = ax.plot(ell, Doppler * scale_3    , ls='solid', color='black' , label="Doppler")



    leg1 = plt.legend(handles=[Ctot], loc='upper right')
    leg2 = plt.legend(handles=[C1, C2, C3, C4], loc='upper left')

    ax.add_artist(leg1)
    ax.add_artist(leg2)

    set_ax_info(ax, xlabel, ylabel, legend=False, 
                ylim=ylim, xlim=xlim, 
                yticks=yticks,
                ypad=ypad)
    if logy:
        ax.set_yscale('log')

    if logx:
        ax.set_xscale('log')
    
    if xticks is not None:
        ticklabels=[f"${tick}$" for tick in xticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(ticklabels)

    
    save_push(fig, fname, save, push, temp=temp)

def plot_matter_PS(k, Pk, k_eq,
                   galaxy_data, wmap_data,
                   fname,  
                   logx=True, logy=True, 
                   xlabel=None, ylabel=None,  
                   xlim=None, ylim=None,  
                   ypad=None, 
                   figsize=(10,6),  
                   save=True,push=False, temp=False):


    fig, ax = plt.subplots(figsize=figsize)

    k_galaxy, Pk_galaxy, error_galaxy = galaxy_data
    k_wmap, Pk_wmap, error_wmap = wmap_data
    # error_wmap = np.abs(error_wmap - Pk_wmap)
    
    P1, = ax.plot(k, Pk, color='blue', label='Prediction')
    P2 = ax.errorbar(k_galaxy, Pk_galaxy, error_galaxy, barsabove=True, fmt='o',
                capthick=1.5, capsize=5, elinewidth=2, color='red', ms=3, 
                label='SDSS Galaxies (DR7 LRG)')
    P3 = ax.errorbar(k_wmap, Pk_wmap, error_wmap, barsabove=True, fmt='o',
                capthick=1.5, capsize=5, elinewidth=2, color='green', ms=3, 
                label='CMB (WMAP+ACT)')

    if ylim is None:
        ylim_vline = ax.get_ylim()
    else:
        ylim_vline = ylim

    if k_eq is not None:
        vl = ax.vlines(k_eq.value, *ylim_vline, colors='black', ls='dashed', label=r"$k_\mathrm{eq}$")

    leg1 = plt.legend(handles=[P1, P2, P3], loc='lower left')
    leg2 = plt.legend(handles=[vl], loc='lower right')

    ax.add_artist(leg1)
    ax.add_artist(leg2)


    
    set_ax_info(ax, xlabel, ylabel, legend=False, 
                ylim=ylim, xlim=xlim, ypad=ypad, xpad=10)

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

    for idx, l in enumerate(ell_values):
        ax.plot(x, y[idx], lw=1, label=fr"$\ell=${str(l)}")




    
    set_ax_info(ax, xlabel, ylabel, legend=True, 
                ylim=ylim, xlim=xlim, yticks=yticks,
                ypad=ypad)
    if logy:
        ax.set_yscale('log')

    if logx:
        ax.set_xscale('log')



    save_push(fig, fname, save, push, temp=temp)


def subplot_Theta0(x, y,
               ell_values,
               x_dec,
               fname,
               logx=True, 
               logy=True,
               color='blue',
               xlabel=r"$x$", 
               ylabel=None, 
               xlim=None, 
               ylim=None, 
               legendloc='best', 
               yticks=None, 
               ypad=None,
               figsize=(10,7), 
               save=False,push=False, temp=False):


    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'hspace': 0})


    if ylim is None:
        y0 = np.min(y) - 0.5
        y1 = np.max(y) + 0.5
        ylim = [y0, y1]
        ax[0].set_ylim(ylim)
        ax[1].set_ylim(ylim)


    x_dec_fill_low = x_dec-0.1
    x_dec_fill_high = x_dec+0.1

    fill = ax[0].fill_betweenx(ylim, x_dec_fill_low, x_dec_fill_high, color='gray', alpha=0.5,
                            label=r"$x\sim x_\mathrm{dec}$")

    fill = ax[1].fill_betweenx(ylim, x_dec_fill_low, x_dec_fill_high, color='gray', alpha=0.5,
                            label=r"$x\sim x_\mathrm{dec}$")

    lines = []
    # lines2 = []
    ls = ['--', '-.', ':', '--', '-.', ':']
    cs = ["red", "red", "red", "blue", "blue", "blue"]
    aaxx = [0,0,0,1,1,1]
    
    for idx, Theta0 in enumerate(y):
        label = rf"$k = {{{ell_values[idx]:.0f}}}/\eta_0 $"
        line, = ax[aaxx[idx]].plot(x, Theta0, color=cs[idx], ls=ls[idx],  label=label)
        # line2, = ax.plot(x, Theta0, color=color, ls=ls[idx],  label=label)

        lines.append(line)

    # leg1 = ax[0].legend(handles=[fill], loc="lower right")
    leg2 = ax[0].legend(handles=lines[0:3], loc="lower left", ncol=3)
    # ax[0].add_artist(leg1)
    # ax.add_artist(leg1)
    ax[0].add_artist(leg2)

    leg1 = ax[1].legend(handles=[fill], loc="upper right")
    leg3 = ax[1].legend(handles=lines[3:], loc="lower left", ncol=3)
    ax[1].add_artist(leg1)


    ax[1].add_artist(leg3)


    
    set_ax_info(ax[0], xlabel, ylabel=None, legend=False, 
                ylim=ylim, xlim=xlim, yticks=yticks,
                ypad=ypad)


    set_ax_info(ax[1], xlabel, ylabel=None, legend=False, 
                ylim=ylim, xlim=xlim, yticks=yticks,
                ypad=ypad)

    plt.subplots_adjust(hspace=0)

    yticks = ax[1].get_yticklabels()
    yticks[-1].set_visible(False)
    ax[0].set_title(ylabel)
    # fig.text(0.02, 0.55, 
            #  ylabel, 
            #  va='center', 
            #  ha='center', 
            #  rotation='vertical', 
            #  fontsize=plt.rcParams['axes.labelsize'])

    save_push(fig, fname, save, push, temp=temp)


def plot_Theta0(x, y,
               ell_values,
               x_dec,
               fname,
               logx=True, 
               logy=True,
               color='blue',
               xlabel=r"$x$", 
               ylabel=None, 
               xlim=None, 
               ylim=None, 
               legendloc='best', 
               yticks=None, 
               ypad=None,
               figsize=(10,7), 
               save=False,push=False, temp=False):


    fig, ax = plt.subplots(figsize=figsize)


    if ylim is None:
        y0 = np.min(y) - 0.5
        y1 = np.max(y) + 0.5
        ylim = [y0, y1]
        ax[0].set_ylim(ylim)
        ax[1].set_ylim(ylim)


    x_dec_fill_low = x_dec-0.1
    x_dec_fill_high = x_dec+0.1

    fill = ax.fill_betweenx(ylim, x_dec_fill_low, x_dec_fill_high, color='gray', alpha=0.5,
                            label=r"$x\sim x_\mathrm{dec}$")

    lines = []
    ls = ['--', '-.', ':']
    
    for idx, Theta0 in enumerate(y):
        label = rf"$k = {{{ell_values[idx]:.0f}}}/\eta_0 $"
        line, = ax.plot(x, Theta0, color=color, ls=ls[idx],  label=label)

        lines.append(line)

    leg1 = ax.legend(handles=[fill], loc="lower right")
    leg2 = ax.legend(handles=lines, loc="lower left", ncol=3)

    ax.add_artist(leg1)
    ax.add_artist(leg2)


    
    set_ax_info(ax, xlabel, ylabel=ylabel, legend=False, 
                ylim=ylim, xlim=xlim, yticks=yticks,
                ypad=ypad)


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