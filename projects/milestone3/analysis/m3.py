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
global RECOMBINATION
SAVE        = False 
TEMP        = False 
RECOMBINATION  = False 


class Perturbations:
    def __init__(self, 
                 filename,
                 f2, f3,
                 length_unit=u.Mpc,
                 time_unit=u.yr):
        
        files = [filename, f2, f3]
        self.k_vals = []
        file_dict = {}
        for file in files:
            # file = file.split(".")
            # print(file)
            if file.startswith("perturbations_k"):
                k_val = float(file.split("_")[1][:-4][1:])
                self.k_vals.append(k_val)
                file_dict[k_val] = file 
        mpc_label = r"\mathrm{Mpc}"
        self.k_labels = [rf"$k={k}/{mpc_label}$" for k in self.k_vals]

        self.data = self.load_data(file_dict)


        self.x              = self.data[0][0]

        self.length_unit    = length_unit
        self.time_unit      = time_unit



    
    def load_data_from_txt(self, filename, skiprows=1):
        return np.loadtxt(data_path + filename, unpack=True, skiprows=skiprows) 

    def load_data(self, file_dict, skiprows=1):
        data = np.array([self.load_data_from_txt(file_dict[k]) for k in self.k_vals])
        data = np.transpose(data, (1,0,2))
        return data 

    
    
    def load_delta(self):
        self.deltas    = True 
        # self.delta_cdm = self.data[1] 
        # self.delta_cdm2 = self.d2[1] 
        # self.delta_cdm3 = self.d3[1] 
        # self.delta_b   = self.data[2]
        # self.delta_b2   = self.d2[2]
        # self.delta_b3   = self.d3[2]
        self.delta_cdm = self.data[1]
        self.delta_b   = self.data[2]
    

    def load_v(self):
        self.vels  = True 
        # self.v_cdm = self.data[3]
        # self.v_b   = self.data[4]
        # self.v_cdm2 = self.d2[3]
        # self.v_b2   = self.d2[4]
        # self.v_cdm3 = self.d3[3]
        # self.v_b3   = self.d3[4]
        self.v_cdm = self.data[3]
        self.v_b   = self.data[4]


    def load_Thetas(self):
        self.Thetas = True 

        # self.Theta0 = self.data[5]
        # self.Theta1 = self.data[6]
        # self.Theta2 = self.data[7]
        # self.Theta0_2 = self.d2[5]
        # self.Theta1_2 = self.d2[6]
        # self.Theta2_2 = self.d2[7]
        # self.Theta0_3 = self.d3[5]
        # self.Theta1_3 = self.d3[6]
        # self.Theta2_3 = self.d3[7]
        self.Theta0 = self.data[5]
        self.Theta1 = self.data[6]
        self.Theta2 = self.data[7]

        
    def load_fields(self):
        self.fields = True 
        self.Phi    = self.data[8]
        self.Psi    = self.data[9]
        self.Pi     = self.data[10]

        # self.Phi_2    = self.d2[8]
        # self.Psi_2    = self.d2[9]
        # self.Pi_2     = self.d2[10]
        # self.Phi_3    = self.d3[8]
        # self.Psi_3    = self.d3[9]
        # self.Pi_3     = self.d3[10]


    def x_to_redshift(self, x):
        # Convert from x to redshift 
        return np.exp(-x) - 1 
    
    ##########################################
    # Make use of the times data in the class 
    ##########################################
    def assert_valid_recombination_value(self, print_x=False, stop=False):
        """
        Check that recombination occurs at z in [1050, 1150]
        both via tau(x_rec)=1 and g_tilde(x_rec)=max(g_tilde)
        """
        if hasattr(self, 'tau'):
            pass
        else:
            self.load_taus()

        if hasattr(self, 'g_tilde'):
            pass
        else:
            self.load_g()

        tau_equal_1_idx = np.argmin(np.abs(self.tau - 1))
        g_peak_idx      = np.argmax(self.g_tilde)
        x_rec_tau = self.x[tau_equal_1_idx]
        x_rec_g= self.x[g_peak_idx]
        z_rec_tau = self.x_to_redshift(x_rec_tau)
        z_rec_g   = self.x_to_redshift(x_rec_g)
    

        if z_rec_tau < 1050 or z_rec_tau > 1150:
            print(" ")
            print("### ERROR ###")
            print("Invlaid solution. tau=1 at z outside [1050, 1150]")
            print_x = True 
            stop = True  

        elif z_rec_g < 1050 or z_rec_g > 1150:
            print("### ERROR ###")
            print("Invlaid solution. g peaks at z outside [1050, 1150]")
            print_x = True 
            stop = True 

        if print_x:
            print(f"z(tau=1) = {z_rec_tau:.5f}, x={x_rec_tau:.5f}")
            print(f"z(g_max) = {z_rec_g:.5f}, x={x_rec_g:.5f}")

        if stop:
            exit()

    def assert_normalized_g_tilde(self, tol=1e-4, print_g=False, stop=False):
        """
        Check that g_tilde is sufficiently normalized to 1
        by checking whether abs(g-1) > tol   
        """

        if hasattr(self, 'g_tilde'):
            pass
        else:
            self.load_g()

        g_integrated = simpson(self.g_tilde, self.x)
        delta_g_norm = abs(g_integrated - 1)
        if delta_g_norm > tol:
            print("### WARNING ###")
            print("g_tilde(x) not properly normalized.")
            print_g = True 
            stop = True
        
        if print_g:
            print(f"  int g(x) deviation = {delta_g_norm:.5e}")

        if stop:
            exit() 

    
    def plot_delta(self, xlim=[-12,0], ylim=[1e-8, 1e6],
                                  figname="tbd.pdf"):

        if hasattr(self, 'deltas'):
            pass
        else:
            self.load_delta()

        # x    = self.x 
        # dcdm, dcdm2, dcdm3 = self.delta_cdm
        # db,db2,db3   = self.delta_b 
        # k1, k2, k3 = self.k_labels
        # db2   = self.delta_b2 
        # dcdm2 = self.delta_cdm2
        # dcdm3 = self.delta_cdm3
        # db3   = self.delta_b3 

        # plt.title(r"$\delta$")
        # plt.plot(x, dcdm,     c='blue', label=k1)
        # plt.plot(x, db, '--', c='blue')
        # plt.plot(x, dcdm2,     c='orange', label=k2)
        # plt.plot(x, db2, '--', c='orange')
        # plt.plot(x, dcdm3,     c='green', label=k3)
        # plt.plot(x, db3, '--', c='green')
        # plt.ylim([1e-1, 1e5])
        # plt.yscale('log')
        # plt.legend()
        # plt.show()
        # ylabel = r"$\delta_\mathrm{CDM},\,\delta_{b}$"
        ylegends = [r"$\delta_\mathrm{CDM}$", r"$\delta_b$"]
        plot.plot_cdm_baryon_for_n_k_values(self.x, self.delta_cdm, self.delta_b, 
                                            self.k_labels,
                                            ylegends, 
                                            fname="deltas.pdf",
                                            ylabel=None, ylim=[1e-1, 1e5],
                                            legendloc1='upper left', 
                                            legendloc2='center left',
                                            save=SAVE, temp=TEMP)


    def plot_v(self, xlim=[-12,0], ylim=[1e-8, 1e6],
                                  figname="tbd.pdf"):

        if hasattr(self, 'vels'):
            pass
        else:
            self.load_v()

        ylegends = [r"$v_\mathrm{CDM}$", r"$v_b$"]
        plot.plot_cdm_baryon_for_n_k_values(self.x, self.v_cdm, self.v_b, 
                                            self.k_labels,
                                            ylegends, 
                                            fname="vels.pdf",
                                            ylabel=None, ylim=[1e-6, 1e2],
                                            legendloc1='upper left', 
                                            legendloc2='lower right',
                                            save=SAVE, temp=TEMP)
        

    
    def plot_Theta(self, xlim=[-12,0], ylim=[1e-8, 1e6],
                                  figname="tbd.pdf"):

        if hasattr(self, 'Thetas'):
            pass
        else:
            self.load_Thetas()

        x      = self.x 
        Theta0 = self.Theta0
        Theta1 = self.Theta1
        Theta2 = self.Theta2

        Theta0_2 = self.Theta0_2
        Theta1_2 = self.Theta1_2
        Theta2_2 = self.Theta2_2
        Theta0_3 = self.Theta0_3
        Theta1_3 = self.Theta1_3
        Theta2_3 = self.Theta2_3
        plt.title(r"$\Theta_0$")
        plt.plot(x, Theta0,   c='blue', label='k=0.001/Mpc')
        plt.plot(x, Theta0_2, c='orange', label='k=0.01/Mpc')
        plt.plot(x, Theta0_3, c='green', label='k=0.1/Mpc')
        plt.legend()
        plt.ylim([-1, 1])
        # plt.yscale('log')
        plt.show()

        plt.title(r"$\Theta_1$")
        plt.plot(x, Theta1,   c='blue', label='k=0.001/Mpc')
        plt.plot(x, Theta1_2, c='orange', label='k=0.01/Mpc')
        plt.plot(x, Theta1_3, c='green', label='k=0.1/Mpc')
        plt.legend()
        # plt.yscale('log')
        plt.ylim([-1, 1])
        plt.show()


    def make_table(self, saha=False):
        if saha:
            fname = self.rec_dec_saha
        else:
            fname = self.rec_dec_file

        x, z, t, r = np.loadtxt(data_path + fname, skiprows=1, usecols=(1,2)) 

        
        t       = (t * u.s).to(self.time_unit).value 
        r       = (r * u.m).to(self.length_unit).value


        plot.time_table(x, z, t, r,
                        saha=saha,
                        save=SAVE,
                        temp=TEMP)








p = Perturbations(filename="perturbations_k0.001.txt",
                  f2="perturbations_k0.01.txt",
                  f3="perturbations_k0.1.txt")
p.plot_delta()
# p.plot_v()
# p.plot_Theta()

