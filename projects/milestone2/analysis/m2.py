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

global SAVE 
global TEMP
global DECOUPLING
SAVE        = False 
TEMP        = False 
DECOUPLING  = False 


class Recombination:
    def __init__(self, 
                 filename,
                 rec_and_dec_times_fname, 
                 saha=False,
                 length_unit=u.Mpc,
                 time_unit=u.Myr):
        
        self.data           = self.load_data(filename)
        self.x              = self.data[0]

        self.rec_dec_file   = rec_and_dec_times_fname

        self.length_unit    = length_unit
        self.time_unit      = time_unit

        self.saha = saha 


    
    def load_data(self, filename, skiprows=2):
        return np.loadtxt(data_path + filename, unpack=True, skiprows=skiprows)
    

    def load_x_decoupling_tau(self):
        self.xdec_tau = np.loadtxt(data_path + self.rec_dec_file, skiprows=1).T[0][0] 

    


    def load_Xe(self):
        self.Xe             = self.data[1]
        

    def load_ne(self, convert_unit=True):
        self.ne     = self.data[2] * u.m**(-3) 
        if convert_unit:
            self.ne = self.ne.to(self.length_unit**(-3))


    def load_taus(self):
        self.tau            = self.data[3]
        self.dtau_dx        = self.data[4]
        self.ddtau_ddx      = self.data[5]


    def load_g(self):
        self.g_tilde        = self.data[6]
        self.dg_tilde_dx    = self.data[7]
        self.ddg_tilde_ddx  = self.data[8]

    def load_sound_horizon(self, convert_unit=True):
        self.sound_horizon = self.data[9] * u.m 
        if convert_unit:
            self.sound_horizon = self.sound_horizon.to(self.length_unit)

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


    def compare_Xe(self, x_saha, Xe_saha, 
                   xdec_peebles=None, xdec_saha=None,
                   ylim=[1e-4, 2], xlim=[-8,-5]):
        
        if hasattr(self, 'Xe'):
            pass
        else:
            self.load_Xe()

        if DECOUPLING:
            if not hasattr(self, "xdec_tau"):
                self.load_x_decoupling_tau()
            
            xdec_peebles = self.xdec_tau

            

        plot.compare_Xe_peebles_and_saha(self.x, x_saha,
                                         self.Xe, Xe_saha,
                                         xdec_peebles=xdec_peebles, xdec_saha=xdec_saha,
                                         fname="compare_Xe_peebles_saha.pdf",
                                         xlim=xlim, ylim=ylim,
                                         decoupling_times=DECOUPLING,
                                         save=SAVE, temp=TEMP)



    def plot_visibility_functions(self, dg_dx_scaling=10, 
                                        ddg_ddx_scaling=300,
                                        xlim=[-7.5, -6],
                                        ylim=None):
        
        """
        Plot g(x), g'(x) and g''(x)
        scale g'(x) and g''(x) to fit in the same plot as g(x)
        """

        if hasattr(self, 'g_tilde'):
            pass
        else:
            self.load_g()

        dg_dx_scaled = self.dg_tilde_dx / dg_dx_scaling
        ddg_ddx_scaled = self.ddg_tilde_ddx / ddg_ddx_scaling

        x   = self.x 
        y   = self.g_tilde
        dy  = dg_dx_scaled 
        ddy = ddg_ddx_scaled

        dg_str = str(dg_dx_scaling)

        y_legend    = r"$\tilde{g}(x)$"
        dy_legend   = r"$\tilde{g}'(x)$" + rf"$/{str(dg_dx_scaling)}$"
        ddy_legend  = r"$\tilde{g}''(x)$" + rf"$/{str(ddg_ddx_scaling)}$"


        plot.plot_quantity_with_derivatives(x, y, dy, ddy, 
                                            y_legend, dy_legend, ddy_legend,
                                            fname="g_plot.pdf",
                                            xlim=xlim, log=False,
                                            save=SAVE, temp=TEMP)


    def plot_tau_with_derivatives(self, xlim=[-10,0], ylim=[1e-8, 1e6]):

        if hasattr(self, 'tau'):
            pass
        else:
            self.load_taus()

        x   = self.x 
        y   = self.tau 
        dy  = - self.dtau_dx
        ddy = self.ddtau_ddx

        y_legend   = r"$\tau(x)$"
        dy_legend  = r"$- \tau'(x)$"
        ddy_legend = r"$\tau''(x)$"

        plot.plot_quantity_with_derivatives(x, y, dy, ddy, 
                                            y_legend, dy_legend, ddy_legend,
                                            fname="tau_plot.pdf",
                                            xlim=xlim, ylim=ylim,
                                            save=SAVE, temp=TEMP)


    def make_table(self):

        rec_dec_times = np.loadtxt(data_path + self.rec_dec_file, skiprows=1).T 

        self.x_times_data  = rec_dec_times[0]
        self.z_times_data  = rec_dec_times[1]
        self.t_times_data  = rec_dec_times[2] * u.s 
        self.rs_times_data = rec_dec_times[3] * u.m 

        self.t_times_data  = self.t_times_data.to(self.time_unit)
        self.rs_times_data = self.rs_times_data.to(self.length_unit)


        plot.time_table(self.x_times_data,
                        self.z_times_data,
                        self.t_times_data,
                        saha=self.saha,
                        save=SAVE,
                        temp=TEMP)






rec = Recombination(filename="recombination.txt", 
                    rec_and_dec_times_fname="rec_times.txt",
                    saha=False)

rec_saha_only = Recombination(filename="recombination_saha.txt", 
                              rec_and_dec_times_fname="rec_times_saha.txt",
                              saha=True)
rec_saha_only.load_Xe()
x_saha, Xe_saha = rec_saha_only.x, rec_saha_only.Xe

rec_saha_only.load_x_decoupling_tau()
x_dec_saha = rec_saha_only.xdec_tau


rec.assert_valid_recombination_value()
rec.assert_normalized_g_tilde()



# SAVE        = True
# TEMP        = True
# DECOUPLING  = True 

# rec.plot_tau_with_derivatives()
# rec.plot_visibility_functions()
# rec.compare_Xe(x_saha=x_saha, Xe_saha=Xe_saha, xdec_saha=x_dec_saha)

# rec.make_table()
# rec_saha_only.make_table()




