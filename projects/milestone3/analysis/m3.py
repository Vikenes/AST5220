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
                 length_unit=u.Mpc,
                 time_unit=u.yr):
        
        self.data           = self.load_data(filename)
        self.x              = self.data[0]


        self.length_unit    = length_unit
        self.time_unit      = time_unit



    
    def load_data(self, filename, skiprows=1):
        return np.loadtxt(data_path + filename, unpack=True, skiprows=skiprows)
    
    def load_delta(self):
        self.delta     = True 
        self.delta_cdm = self.data[1] 
        self.delta_b   = self.data[2]
    
    def load_v(self):
        self.v     = True 
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

        if hasattr(self, 'delta'):
            pass
        else:
            self.load_delta()

        x    = self.x 
        dcdm = self.delta_cdm
        db   = self.delta_b 

        plt.plot(x, dcdm, label='cdm')
        plt.plot(x, db, '--', label='b')
        plt.yscale('log')
        plt.show()


    def plot_v(self, xlim=[-12,0], ylim=[1e-8, 1e6],
                                  figname="tbd.pdf"):

        if hasattr(self, 'v'):
            pass
        else:
            self.load_v()

        x    = self.x 
        vcdm = self.v_cdm
        vb   = self.v_b 

        plt.plot(x, vcdm, label='cdm')
        plt.plot(x, vb, '--', label='b')
        plt.yscale('log')
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






p = Perturbations(filename="perturbations_k0.001.txt")
p.plot_v()
