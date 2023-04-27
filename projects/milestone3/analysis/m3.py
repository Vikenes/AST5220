import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
import os 
import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning )

from scipy.integrate import simpson 


"""
All calculations and such are performed in this file,
but the actual plotting is done in plot.py
"""
import plot

data_path = "/home/vetle/Documents/master_studies/subjects/V23/AST5220/projects/milestone3/data/"

global SAVE 
global TEMP
global XLIM
SAVE        = False 
TEMP        = False 
XLIM        = [-15, 0]

class Perturbations:
    def __init__(self, 
                 f1, f2, f3,
                 length_unit=u.Mpc,
                 time_unit=u.yr):
        
        ### Ensure correct k-value for each file 
        ### by making dictionary with keys from k-values in filename 
        files           = [f1, f2, f3]
        self.file_dict  = self.make_file_dictionary(files) 
        
        # Make labels for plotting 
        mpc_label     = r"\mathrm{Mpc}"
        self.k_labels = [rf"$k={k}/{mpc_label}$" for k in self.k_vals]


        self.data           = self.load_data(self.file_dict)
        self.x              = self.data[0][0]

        self.length_unit    = length_unit
        self.time_unit      = time_unit

        self.tc_end = np.array([float(np.loadtxt(data_path + f, max_rows=1)) for f in files])

        # self.tc_end = np.array(list(set(tc_end)))


    def make_file_dictionary(self, filenames):
        self.k_vals = []
        file_dict = {}
        for file in filenames:
            if file.startswith("perturbations_k"):
                k_val = float(file.split("_")[1][:-4][1:])
                self.k_vals.append(k_val)
                file_dict[k_val] = file 

        return file_dict 
    
    def load_data_from_txt(self, filename, skiprows=1):
        return np.loadtxt(data_path + filename, unpack=True, skiprows=skiprows) 

    def load_data(self, file_dict, skiprows=1):
        data = np.array([self.load_data_from_txt(file_dict[k]) for k in self.k_vals])
        data = np.transpose(data, (1,0,2))

        return data 

    
    
    def load_delta(self):
        self.deltas    = True 
        
        self.delta_cdm = self.data[1]
        self.delta_b   = self.data[2]

        # print(np.min(self.delta_cdm,axis=1))
        # print(np.min(self.delta_b,axis=1))

    

    def load_v(self):
        self.vels  = True 
        
        self.v_cdm = self.data[3]
        self.v_b   = self.data[4]


    def load_Thetas(self):
        self.Thetas = True 

        self.Theta0 = self.data[5]
        self.Theta1 = self.data[6]
        self.Theta2 = self.data[7]

        
    def load_fields(self):
        self.fields = True 
        self.Phi    = self.data[8]
        self.Psi    = self.data[9]
        self.Pi     = self.data[10]
    
    def load_source_func(self):
        self.source         = True 
        self.S_tilde        = self.data[11]
        self.S_tilde_j5     = self.data[12]
        self.S_tilde_j50    = self.data[13] 
        self.S_tilde_j500   = self.data[14] 


    def plot_delta(self, xlim=XLIM, ylim=[1e-2, 1e4]):

        if hasattr(self, 'deltas'):
            pass
        else:
            self.load_delta()
        ylabel = r"$\delta_\mathrm{CDM},\,\delta_b$"
        ylegends = ["CDM", "baryons"]
        plot.plot_cdm_baryon_for_n_k_values(self.x, self.delta_cdm, self.delta_b, 
                                            self.k_labels,
                                            ylegends, 
                                            fname="deltas",
                                            xend=self.tc_end,
                                            xlim=xlim,
                                            ylabel=ylabel, ylim=ylim,
                                            legendloc1='upper left', 
                                            legendloc2='lower right',
                                            legendloc3='lower left',
                                            figsize=(10,6),
                                            save=SAVE, temp=TEMP)

    def plot_delta_gamma(self, xlim=XLIM, ylim=[-2, 4],
                         log=False, yabs=False):

        if hasattr(self, 'Thetas'):
            pass
        else:
            self.load_Thetas()

        delta_gamma = 4 * self.Theta0
        ylabel = r"$\delta_\gamma$"

        plot.plot_quantity_for_n_k_values(x=self.x, y=delta_gamma, 
                                          k_legends = self.k_labels, 
                                          fname="delta_gamma",
                                          xend=self.tc_end,
                                          absval_of_y=yabs,
                                          xlabel=r"$x$", ylabel=ylabel,
                                          ylim=ylim, log=log,
                                          xlim=xlim,
                                          save=SAVE, temp=TEMP)


    def plot_v(self, xlim=XLIM, ylim=[1e-5, 1e3]):

        if hasattr(self, 'vels'):
            pass
        else:
            self.load_v()

        ylabel = r"$v_\mathrm{CDM},\, v_b$"
        ylegends = ["CDM", "baryons"]
        plot.plot_cdm_baryon_for_n_k_values(self.x, self.v_cdm, self.v_b, 
                                            self.k_labels,
                                            ylegends, 
                                            fname="vels",
                                            xend=self.tc_end,
                                            ylabel=ylabel, ylim=ylim,
                                            xlim=xlim,
                                            legendloc1='upper left', 
                                            legendloc2='lower right',
                                            legendloc3='upper right',
                                            figsize=(10,6),
                                            save=SAVE, temp=TEMP)
        
    def plot_v_gamma(self, xlim=XLIM, ylim=[-2, 4], log=False, yabs=False):

        if hasattr(self, 'Thetas'):
            pass
        else:
            self.load_Thetas()

        v_gamma = -3 * self.Theta1
        ylabel = r"$v_\gamma$"

        plot.plot_quantity_for_n_k_values(x=self.x, y=v_gamma, 
                                          k_legends = self.k_labels, 
                                          fname="v_gamma",
                                          xend=self.tc_end,
                                          absval_of_y=yabs,
                                          xlabel=r"$x$", ylabel=ylabel,
                                          ylim=ylim, log=log,
                                          xlim=xlim,
                                          legendloc='upper left',
                                          save=SAVE, temp=TEMP)
    
    def plot_Theta(self, Theta_number=0, xlim=XLIM, ylim=[-0.5, 1], legendloc='best'):

        if hasattr(self, 'Thetas'):
            pass
        else:
            self.load_Thetas()

        thetas = [self.Theta0, self.Theta1, self.Theta2]
        theta = thetas[Theta_number]
        theta_label = fr"$\Theta_{Theta_number}$"
        theta_fname = f"Theta{Theta_number}"

        plot.plot_quantity_for_n_k_values(x=self.x, y=theta, 
                                          k_legends = self.k_labels, 
                                          fname=theta_fname,
                                          xend=self.tc_end,
                                          xlabel=r"$x$", ylabel=theta_label,
                                          ylim=ylim, log=False,
                                          xlim=xlim,
                                          legendloc=legendloc,
                                          save=SAVE, temp=TEMP)

    def plot_Phi(self, xlim=XLIM, ylim=[0, 0.7]):

        if hasattr(self, 'fields'):
            pass
        else:
            self.load_fields()


        plot.plot_quantity_for_n_k_values(x=self.x, y=self.Phi, 
                                          k_legends = self.k_labels, 
                                          fname="Phi",
                                          xend = self.tc_end,
                                          xlabel=r"$x$", ylabel=r"$\Phi$",
                                          ylim=ylim, xlim=xlim, 
                                          log=False,
                                          save=SAVE, temp=TEMP)
    
    def plot_Phi_plus_Psi(self, xlim=XLIM, ylim=None):

        if hasattr(self, 'fields'):
            pass
        else:
            self.load_fields()

        y = self.Phi + self.Psi
        ylabel = r"$\Phi + \Psi$"
        plot.plot_quantity_for_n_k_values(x=self.x, y=y, 
                                          k_legends = self.k_labels, 
                                          fname="Phi_plus_Psi",
                                          xend=self.tc_end,
                                          xlabel=r"$x$", ylabel=ylabel,
                                          ylim=ylim, xlim=xlim, 
                                          log=False,
                                          legendloc='upper left',
                                          save=SAVE, temp=TEMP)

    def plot_source_function(self, xlim=XLIM, ylim=None, no=0, fname=None):
        if hasattr(self, 'source'):
            pass
        else:
            self.load_source_func()

        if no == 0:
            if fname is None:
                fname = "S_alone"
            plot.plot_source_function(self.x, self.S_tilde,
                                      k_legends=self.k_labels,
                                      ylabel=r"$\tilde{S}$",
                                      fname=fname,
                                      xlim=xlim,ylim=ylim,
                                      save=SAVE, temp=TEMP)

        if no == 1:
            if fname is None:
                fname = "S_j5"
            plot.plot_source_function(self.x, self.S_tilde_j5, 
                                      k_legends=self.k_labels,
                                      ylabel=r"$\tilde{S} j_5 $",
                                      fname=fname,
                                      xlim=xlim,ylim=ylim,
                                      save=SAVE, temp=TEMP)

        elif no == 2:
            if fname is None:
                fname = "S_j50"
            plot.plot_source_function(self.x, self.S_tilde_j50, 
                                      k_legends=self.k_labels,
                                      ylabel=r"$\tilde{S} j_{50} $",
                                      fname=fname,
                                      xlim=xlim,ylim=ylim,
                                      save=SAVE, temp=TEMP)
            # plt.title(r"$\tilde{S}(k,x) j_{50}$")
            # plt.plot(self.x, self.S_tilde_j50[2])
            # plt.xlim(xlim)
            # plt.ylim(ylim)
            # plt.show()
        elif no == 3:
            if fname is None:
                fname = "S_j500"
            plot.plot_source_function(self.x, self.S_tilde_j500, 
                                      k_legends=self.k_labels,
                                      ylabel=r"$\tilde{S} j_{500} $",
                                      fname=fname,
                                      xlim=xlim,ylim=ylim,
                                      save=SAVE, temp=TEMP)
            # plt.title(r"$\tilde{S}(k,x) j_{500}$")
            # plt.plot(self.x, self.S_tilde_j500[2])
            # plt.xlim(xlim)
            # plt.ylim(ylim)
            # plt.show()
        else:
            pass     

   





p = Perturbations(f1="perturbations_k0.001.txt",
                  f2="perturbations_k0.01.txt",
                  f3="perturbations_k0.1.txt")

SAVE=True
# TEMP=True
# p.plot_delta()
# p.plot_delta_gamma(ylim=[-3, 5])
# p.plot_v()
# p.plot_v_gamma(ylim=[-1.5,1.5])

# p.plot_Theta(0, ylim=[-0.8,1], legendloc='lower left')
# p.plot_Theta(1, ylim=[-0.5, 0.6], legendloc='upper left')
p.plot_Theta(2, xlim=[-10,0], ylim=[-0.1, 0.2], legendloc='upper right')

# p.plot_Phi()
# p.plot_Phi_plus_Psi(ylim=[-0.006,0.026])

#### NOT USING ###### 
# p.plot_delta_gamma(ylim=[1e-2,1e1], yabs=True, log=True) ## Not using 

# p.plot_source_function(xlim=[-7.5,-6], no=1, fname='S_j5_center.pdf')
# p.plot_source_function(xlim=[-1,0], ylim=[-0.026,0.055], no=1, fname='S_j5_end.pdf')

# p.plot_source_function(xlim=[-7.5,-6], no=2, fname='S_j50_center.pdf')
# p.plot_source_function(xlim=[-3,0], ylim=[-0.001,0.0015], no=2, fname='S_j50_mid_end.pdf')
# p.plot_source_function(xlim=[-1,0], ylim=[-0.001,0.0015], no=2, fname='S_j50_end.pdf')

# p.plot_source_function(xlim=[-7.5,-6], no=3, fname='S_j500_center.pdf')
# p.plot_source_function(xlim=[-5,0], ylim=[-3e-5,6e-5], no=3, fname='S_j500_end.pdf')

# p.plot_source_function(xlim=[-7.5,-6], no=0, fname='S_no_l_center.pdf')
# p.plot_source_function(xlim=[-1,0], ylim=[-0.01,0.5], no=0, fname='S_no_l_end.pdf')
