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


def compare_gs(dg_dx_scaling=30, 
                ddg_ddx_scaling=500,
                xlim=[-8, -6]):
        
        """
        Plot g(x), g'(x) and g''(x)
        scale g'(x) and g''(x) to fit in the same plot as g(x)
        """
        rec1 = Recombination("recombination.txt")
        rec2 = Recombination("rec2.txt")
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


def compare_taus(xlim=[-8,-6]):
        rec1 = Recombination("rec2.txt")
        rec2 = Recombination("recombination.txt")
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






# rec = Recombination("recombination.txt")
rec = Recombination("r.txt")

rec_saha_only = Recombination("recombination_saha.txt")
x_saha, Xe_saha = rec_saha_only.x, rec_saha_only.Xe
rec.assert_valid_solution()
# rec.compare_Xe(x_saha, Xe_saha)

# compare_gs()
# compare_taus()    

# rec_sound.plotfrac()
# rec_sound.plot_tau_with_derivatives(xlim=[-7.5,-6.5])
# rec.plot_visibility_functions()


# z1 = rec_sound.x_to_redshift(rec_sound.x)
# z2 = rec_saha.x_to_redshift(rec_saha.x)

# plt.plot(z1, rec_sound.Xe, label="peb")
# plt.plot(z2, rec_saha.Xe, '--', label="saha")
# plt.xscale('log')
# plt.yscale('log')
# plt.xlim(5e3,1e2)
# plt.legend()
# plt.show()
exit()




save = False # If False, the figures produced are only displayed, not saved. 

data1 = plot.load("rec.txt", skiprows=2)
# data1 = plot.load("rec_rk8.txt", skiprows=2)

# pth = os.path.abspath(".") + "/../data/"
# header1 = np.loadtxt(pth+"recombination_sound.txt", dtype=str, unpack=False, max_rows=1)
# print(header1)
# exit()
# data2 = plot.load("recombination_1e4_1e5.txt", skiprows=2)
# data2 = plot.load("recombination_newode.txt", skiprows=2)

x1 = data1[0]
# xmin = np.abs(x1 + 12).argmin()
# xmax = np.abs(x1).argmin()
# x1 = x1[xmin:xmax]
# data1 = data1[:,xmin:xmax]

# x2 = data2[0]

def load_xe_ne(data):

    Xe = data[1]
    ne = data[2] * u.m**(-3)
    return Xe, ne 

def load_taus(data):
    tau, dtau, ddtau = data[3:6]
    return tau, dtau, ddtau 

def load_gs(data):
    g, dg, ddg = data[6:9]
    return g, dg, ddg 

def load_s(data):
    s = data1[9]
    return s

# g2,dg2,ddg2 = load_gs(data2)

# t2, dt2, ddt2 = load_taus(data2)


def plot_tau():
    t1, dt1, ddt1 = load_taus(data1)
    plt.plot(x1, t1)
    plt.plot(x1, -dt1, '--')
    plt.plot(x1, ddt1, ':')
    plt.xlim(-7.5, -6.5)
    plt.ylim(0,2)
    # plt.yscale('log')
    plt.show()

def plot_g(order=None):
    g1,dg1,ddg1 = load_gs(data1)
    
    if order==0 or order==None:
        plt.plot(x1, g1)
        plt.show()
    if order==1 or order==None:
        plt.plot(x1, dg1)#, '--')
        plt.show()
    if order==2 or order==None:
        plt.plot(x1, ddg1)#, ':')
        plt.show()


g1 = load_gs(data1)[0]
xm8 = np.abs(x1 + 7.5).argmin()
xm6 = np.abs(x1 + 7).argmin()

print(simpson(g1, x1), "\n")

exit()

# plt.plot(x1[xm8:xm6], g1[xm8:xm6])
# plt.show()


def plot_s():
    plt.plot(x1,load_s(data1))
    plt.show()

g = load_gs(data1)[0]
s = load_s(data1)
decoupling = np.argmax(g)
x_decoupling = x1[decoupling]
s_dec = s[decoupling]
z_dec = np.exp(-x_decoupling) - 1

tau = load_taus(data1)[0]
tau1_idx = np.abs(tau - 1).argmin()
z_dec2 = np.exp(-x1[tau1_idx]) - 1 

# print("Vetle: ", z_dec2)
# print("Nanna: ", np.exp(-(-6.9874)) - 1);exit()
print(s_dec)
print(z_dec)
print(z_dec2)

