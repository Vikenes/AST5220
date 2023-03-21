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

# data_path = os.path.abspath(".") + "/../data/"
# data_path = os.path.abspath(".") + "/../data/"
data_path = "/home/vetle/Documents/master_studies/subjects/V23/AST5220/projects/milestone2/data/"




class Recombination:
    def __init__(self, filename, length_units=u.Mpc):
        self.data           = self.load_data(filename)
        self.x              = self.data[0]
        self.Xe             = self.data[1]
        self.ne             = self.data[2] * u.Gpc**(-3)
        self.tau            = self.data[3]
        self.dtau_dx        = self.data[4]
        self.ddtau_ddx      = self.data[5]
        self.g_tilde        = self.data[6]
        self.dg_tilde_dx    = self.data[7]
        self.ddg_tilde_ddx  = self.data[8]
        self.sound_horizon  = self.data[9] * u.m 

        self.ne             = self.ne.to(length_units**(-3))
        self.sound_horizon  = self.sound_horizon.to(length_units)

    
    def load_data(self, filename, skiprows=2):
        return np.loadtxt(data_path + filename, unpack=True, skiprows=skiprows)
    

    def x_to_redshift(self, x):
        return np.exp(-x) - 1 
    
    def assert_valid_solution(self, print_x=False, stop=False):
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




def compare_gs(r1, r2,
                dg_dx_scaling=30, 
                ddg_ddx_scaling=500,
                xlim=[-8, -6]):
        
        rec1 = r1
        rec2 = r2
        x1 = rec1.x 
        x2 = rec2.x 
        g1, dg1, ddg1 = rec1.g_tilde, rec1.dg_tilde_dx/dg_dx_scaling, rec1.ddg_tilde_ddx/ddg_ddx_scaling 
        g2, dg2, ddg2 = rec2.g_tilde, rec2.dg_tilde_dx/dg_dx_scaling, rec2.ddg_tilde_ddx/ddg_ddx_scaling

        plt.plot(x1, g1  , label="old")
        plt.plot(x1, dg1 , ls='solid', label="old")
        plt.plot(x1, ddg1, ls='solid', label="old")
        plt.plot(x2, g2  , ls='solid' ,alpha=0.5,  label="new")
        plt.plot(x2, dg2 , ls='dashed',alpha=0.5, label="new")
        plt.plot(x2, ddg2, ls='dotted',alpha=0.5, label="new")


        plt.legend()
        plt.xlim(xlim)
        plt.show()


def compare_taus(r1, r2, xlim=[-8,-6]):
        rec1 = r1
        rec2 = r2
        x1 = rec1.x 
        x2 = rec2.x 
        t1, dt1, ddt1 = rec1.tau, rec1.dtau_dx, rec1.ddtau_ddx 
        t2, dt2, ddt2 = rec2.tau, rec2.dtau_dx, rec2.ddtau_ddx

        plt.plot(x1, t1  , ls='solid' , label=r"$\tau(x)$")
        plt.plot(x1, -dt1, ls='solid', label=r"$-\tau'(x)$")
        plt.plot(x1, ddt1, ls='solid', label=r"$-\tau'(x)$")
        plt.plot(x2, t2  , ls='solid' ,alpha=0.5, label="old")
        plt.plot(x2, -dt2, ls='dashed',alpha=0.5, label="old")
        plt.plot(x2, ddt2, ls='dotted',alpha=0.5, label="old")
        
        plt.legend()
        # plt.xlim(xlim)
        plt.yscale('log')
        plt.show()






rec1 = Recombination("recombination.txt")
rec2 = Recombination("recombination.txt")

# rec1.assert_valid_solution(True)
# rec2.assert_valid_solution(True)


# compare_gs(rec1, rec2)
# compare_taus(rec1, rec2)    


