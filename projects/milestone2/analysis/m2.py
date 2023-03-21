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
        # Convert from x to redshift 
        return np.exp(-x) - 1 
    
    def assert_valid_recombination_value(self, print_x=False, stop=False):
        """
        Check that recombination occurs at z in [1050, 1150]
        both via tau(x_rec)=1 and g_tilde(x_rec)=max(g_tilde)
        """
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


    def compare_Xe(self, x_saha, Xe_saha):
        
        plt.plot(self.x, self.Xe, label='Saha+peebles')
        plt.plot(x_saha, Xe_saha, ls='dashed', label='saha only')

        plt.ylim(1e-4, 2)
        plt.xlim(-10,0)
        plt.legend()
        plt.yscale('log')
        plt.show()


    def plot_visibility_functions(self, dg_dx_scaling=30, 
                                        ddg_ddx_scaling=500,
                                        xlim=[-7.5, -6.5]):
        
        """
        Plot g(x), g'(x) and g''(x)
        scale g'(x) and g''(x) to fit in the same plot as g(x)
        """

        dg_dx_scaled = self.dg_tilde_dx / dg_dx_scaling
        ddg_ddx_scaled = self.ddg_tilde_ddx / ddg_ddx_scaling

        plt.plot(self.x, self.g_tilde  , label=r"$\tilde{g}(x)$")
        plt.plot(self.x, dg_dx_scaled  , ls='dashed', label=r"$\tilde{g}'(x)$" + f"/{dg_dx_scaling:.0f}")
        plt.plot(self.x, ddg_ddx_scaled, ls='dotted', label=r"$\tilde{g}''(x)$"+ f"/{ddg_ddx_scaling:.0f}")
        plt.legend()
        plt.xlim(xlim)
        plt.show()


    def plot_tau_with_derivatives(self, xlim=[-12,0]):
        
        plt.plot(self.x, self.tau       , ls='solid', label=r"$\tau(x)$")
        plt.plot(self.x, -self.dtau_dx  , ls='dashed', label=r"$-\tau'(x)$")
        plt.plot(self.x, self.ddtau_ddx , ls='dotted', label=r"$\tau''(x)$")
        plt.legend()
        plt.xlim(xlim)
        plt.yscale('log')
        plt.show()



class recomb_and_decoupling_times:
    def __init__(self, filename, length_unit=u.Mpc, time_unit=u.Myr):

        data = np.loadtxt(data_path + filename, skiprows=1).T 
        
        x_values  = data[0]
        z_values  = data[1]
        t_values  = data[2] * u.s 
        rs_values = data[3] * u.m

        t_values  = t_values.to(time_unit)
        rs_values = rs_values.to(length_unit) 


        self.x_dec_tau, self.x_dec_gmax, self.x_rec     = x_values
        self.z_dec_tau, self.z_dec_gmax, self.z_rec     = z_values
        self.t_dec_tau, self.t_dec_gmax, self.t_rec     = t_values
        self.rs_dec_tau, self.rs_dec_gmax, self.rs_rec  = rs_values

        # x_dec_tau, z_dec_tau, t_dec_tau, rs_dec_tau = self.data[0]
        # x_dec_gmax, z_dec_gmax, t_dec_gmax, rs_dec_gmax = self.data[1]
        # x_rec, z_rec, t_rec, rs_rec = self.data[2]
        print(self.x_dec_tau, self.x_dec_gmax, self.x_rec)

        print(self.rs_dec_gmax)
        print(self.t_rec)
    # def load_data_




rec = Recombination("recombination.txt")

rec_saha_only = Recombination("recombination_saha.txt")

# x_saha, Xe_saha = rec_saha_only.x, rec_saha_only.Xe
rec.assert_valid_recombination_value()
rec.assert_normalized_g_tilde()
# rec.plot_visibility_functions()
# rec.compare_Xe(x_saha, Xe_saha)

rec_times = recomb_and_decoupling_times("rec_times.txt")