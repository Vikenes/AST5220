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
global RECOMBINATION
SAVE        = False 
TEMP        = False 
RECOMBINATION  = False 


class Recombination:
    def __init__(self, 
                 filename,
                 rec_and_dec_times_fname=None, 
                 length_unit=u.Mpc,
                 time_unit=u.yr):
        
        self.data           = self.load_data(filename)
        self.x              = self.data[0]

        self.rec_dec_file   = rec_and_dec_times_fname
        self.rec_dec_saha   = rec_and_dec_times_fname.strip(".txt") + "_saha.txt"

        self.length_unit    = length_unit
        self.time_unit      = time_unit

        self.load_taus()
        self.load_g()
        self.load_x_decoupling()
        self.load_Xe()
        self.load_ne()
        self.load_sound_horizon()



    
    def load_data(self, filename, skiprows=2):
        return np.loadtxt(data_path + filename, unpack=True, skiprows=skiprows)
    
    def load_x_decoupling(self):
        x = np.loadtxt(data_path + self.rec_dec_file, skiprows=1, usecols=(1,2,3))[0]
        self.xdec_tau     = x[0]
        self.xdec_g       = x[1]
        self.xrec         = x[2]

    
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
    

    def Xe_today(self):
        self.load_Xe()
        Xe0 = self.Xe[-1]
        print(self.x[-1])
        print(f"X_e at x={self.x[-1]:.2f} is: {self.Xe[-1]:.5e}")

    def compare_Xe(self, x_saha, Xe_saha, 
                   xrec_saha=None,
                   ylim=[1e-4, 2], xlim=[-7.8,-5.3],
                   figname="compare_Xe_peebles_saha.pdf"):
        

        # if RECOMBINATION:
            
            # xrec_peebles = self.xrec


        plot.compare_Xe_peebles_and_saha(self.x, x_saha,
                                         self.Xe, Xe_saha,
                                         xrec_peebles=self.xrec, xrec_saha=xrec_saha,
                                         fname=figname,
                                         xlim=xlim, ylim=ylim,
                                         rec_times=True,
                                         save=SAVE, temp=TEMP)



    def plot_visibility_functions(self, dg_dx_scaling=10, 
                                        ddg_ddx_scaling=300,
                                        xlim=[-7.5, -6],
                                        ylim=None,
                                        figname="g_plot.pdf"):
        
        """
        Plot g(x), g'(x) and g''(x)
        scale g'(x) and g''(x) to fit in the same plot as g(x)
        """

        dg_dx_scaled = self.dg_tilde_dx / dg_dx_scaling
        ddg_ddx_scaled = self.ddg_tilde_ddx / ddg_ddx_scaling

        x   = self.x 
        y   = self.g_tilde
        dy  = dg_dx_scaled 
        ddy = ddg_ddx_scaled
       

        y_legend    = r"$\tilde{g}(x)$"
        dy_legend   = r"$\tilde{g}'(x)$" + rf"$/{str(dg_dx_scaling)}$"
        ddy_legend  = r"$\tilde{g}''(x)$" + rf"$/{str(ddg_ddx_scaling)}$"


        plot.plot_quantity_with_derivatives(x, y, dy, ddy, 
                                            y_legend, dy_legend, ddy_legend,
                                            fname=figname,
                                            x_dec=self.xdec_tau, dec_leg=r"$\tau=1$",
                                            xlim=xlim, ylim=[-3.6, 5.4], log=False,
                                            save=SAVE, temp=TEMP)


    def plot_tau_with_derivatives(self, xlim=[-10,0], ylim=[1e-8, 1e6],
                                  figname="tau_plot.pdf"):

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
                                            x_dec=self.xdec_tau, dec_leg=r"$\tau=1$",
                                            xlim=xlim, ylim=ylim,
                                            save=SAVE, temp=TEMP)


    def make_table(self, saha=False):
        if saha:
            fname = self.rec_dec_saha
        else:
            fname = self.rec_dec_file

        x, z, t, r = np.loadtxt(data_path + fname, skiprows=1, usecols=(1,3)) 

        
        t       = (t * u.s).to(self.time_unit).value 
        r       = (r * u.m).to(self.length_unit).value


        plot.time_table(x, z, t, r,
                        saha=saha,
                        save=SAVE,
                        temp=TEMP)






rec = Recombination(filename="recombination.txt", 
                    rec_and_dec_times_fname="rec_times.txt"
                    )

rec_saha_only = Recombination(filename="recombination_saha.txt", 
                              rec_and_dec_times_fname="rec_times_saha.txt"
                              )
rec_saha_only.load_Xe()
x_saha, Xe_saha = rec_saha_only.x, rec_saha_only.Xe

rec_saha_only.load_x_decoupling()

x_rec_saha = rec_saha_only.xrec




SAVE        = True
# TEMP        = True
# RECOMBINATION  = True 

# rec.Xe_today()

# rec.plot_tau_with_derivatives()
# rec.plot_visibility_functions()
# rec.compare_Xe(x_saha=x_saha, Xe_saha=Xe_saha, xrec_saha=x_rec_saha)
rec.make_table(True)
rec.make_table()

exit()

# rec.plot_tau_with_derivatives(figname="tau_plot_split.pdf")
# rec.plot_visibility_functions(figname="g_plot_split.pdf")
rec.compare_Xe(x_saha=x_saha, 
               Xe_saha=Xe_saha, 
               xdec_saha=x_dec_saha, 
               figname="compare_Xe_peebles_saha_split.pdf")



