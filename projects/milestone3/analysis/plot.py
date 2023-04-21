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
project_path = "/home/vetle/Documents/master_studies/subjects/V23/AST5220/projects/milestone3"

data_path = project_path + "/data/"
latex_path = project_path + "/report/tables/"
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


def set_ax_info(ax, xlabel, ylabel=False, title=None, legend=True, 
                double_legends=[], legendloc='best', legend_size=False):
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




# -----------------------------------------------------------------------------
#   General plot code above
#   XXX
#   Specific plot code below
# -----------------------------------------------------------------------------



def load(file, folder=data_path, skiprows=0):
    # Tired of repeating unpack, delimiter, skiprows for all... 
    return np.loadtxt(folder + file, unpack=True, skiprows=skiprows)



def plot_quantity_for_n_k_values(x, y, k,
                                   k_legends, 
                                   fname, 
                                   xlabel=r"$x$", 
                                   ylabel=None, 
                                   xlim=None, 
                                   ylim=None, 
                                   legendloc='best', 
                                   yticks=None, 
                                   figsize=(8,6), 
                                   log=True, 
                                   save=True, 
                                   temp=False):


    fig, ax = plt.subplots(figsize=figsize)
    y1, y2, y3 = y 
    k1, k2, k3 = k 
    k1_leg, k2_leg, k3_leg = k_legends
    ax.plot(x, y1, color='blue'  , label=k1_leg)
    ax.plot(x, y2, color='orange', label=k2_leg)
    ax.plot(x, y3, color='green' , label=k3_leg)
    
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    if log:
        ax.set_yscale('log')
    
    set_ax_info(ax, xlabel, ylabel, legendloc=legendloc)

    if yticks is not None:
        ax.set_yticks(yticks)

    save_push(fig, fname, save, temp=temp)
    
def plot_cdm_baryon_for_n_k_values(x, y_cdm, y_baryon,
                                   k_legends, ylegends,
                                   fname, 
                                   xlabel=r"$x$", 
                                   ylabel=None, 
                                   xlim=None, 
                                   ylim=None, 
                                   legendloc1='best',
                                   legendloc2='best', 
                                   yticks=None, 
                                   figsize=(8,6), 
                                   log=True, 
                                   save=True, 
                                   temp=False):


    fig, ax = plt.subplots(figsize=figsize)
    y_cdm1   , y_cdm2   , y_cdm3    = y_cdm
    y_baryon1, y_baryon2, y_baryon3 = y_baryon
    y_cdm_legend, y_baryon_legend = ylegends 

    # k1, k2, k3 = k 
    k1_leg, k2_leg, k3_leg = k_legends
    ax.plot(x, y_cdm1   , ls='solid' , color='blue'  , label=k1_leg)
    ax.plot(x, y_baryon1, ls='dashed', color='blue'  )
    ax.plot(x, y_cdm2   , ls='solid' , color='orange', label=k2_leg)
    ax.plot(x, y_baryon2, ls='dashed', color='orange')
    ax.plot(x, y_cdm3   , ls='solid' , color='green' , label=k3_leg)
    ax.plot(x, y_baryon3, ls='dashed', color='green' )


    # ax.plot([], ls='solid' , c='k', label=y_cdm_legend)
    # ax.plot([], ls='dashed', c='k', label=y_baryon_legend)


    
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    if log:
        ax.set_yscale('log')
    set_ax_info(ax, xlabel, ylabel, legendloc=legendloc1)
    # leg1 = ax.legend(loc=legendloc1)
    plt.gca().add_artist(ax.legend(loc=legendloc1))
    solid_line,  = ax.plot([], label=y_cdm_legend   , c='k', linestyle="-")
    dashed_line, = ax.plot([], label=y_baryon_legend, c='k', linestyle="--")
    leg2 = ax.legend([solid_line, dashed_line], ylegends, loc=legendloc2)
    plt.gca().add_artist(leg2)

    if yticks is not None:
        ax.set_yticks(yticks)

    save_push(fig, fname, save, temp=temp)
    




def time_table(x, z, t, r, saha=False, save=False, temp=False):
   

    if saha:
        method = "Saha"
    else:
        method = "Peebles"

    col_labels = [method, r"$x$", r"$z$", r"$t\,\mathrm{[yr]}$"]
    col_fmt = "l|ccc"

    if not saha:
        col_labels.append(r"$r_s \,\mathrm{[Mpc]}$")
        col_fmt += "c"

    row_labels = [r"Decoupling",
                  r"Recombination"]
    
    t1_ = str(np.round(t[0]))
    t2_ = str(np.round(t[1]))

    t1 = t1_[0:3] + "\," + t1_[3:6]
    t2 = t2_[0:3] + "\," + t2_[3:6]

    t = [t1, t2]

    # x = np.round(x,3)
    z = np.round(z,1)

    if saha:
        data = np.array([row_labels, x, z, t]).T 
    else:
        r = [np.round(r[0],2), " "]
        

        data = np.array([row_labels, x, z, t, r]).T 


    table_caption   = "Times when decoupling and recombination occurs, computed from the " 
    table_caption   += method + " equation."

    table_name      = "rec_and_dec_time_table_" + method
    table_label     = "tab:M2:results:" + table_name
    table_fname     = table_name + ".tex"

    if save:
        if temp:
            table_fname = "TEMP" + table_fname

        buffer = latex_path + table_fname 
    else:
        buffer = None 


    table = pd.DataFrame(data).to_latex(
            index=False, 
            header=col_labels, 
            escape=False, 
            column_format=col_fmt,
            caption=table_caption,
            label=table_label,
            position="h",
            buf=buffer
        )


    if save:
        print(f"Saving time table for {method}:")
        print(f"   {table_fname}")
        
    else:
        print(table)


if __name__=='__main__':
    pass 