#ifndef _BACKGROUNDCOSMOLOGY_HEADER
#define _BACKGROUNDCOSMOLOGY_HEADER
#include <iostream>
#include <fstream>
#include "Utils.h"

using Vector = std::vector<double>;

class BackgroundCosmology{
  private:
   
    // Cosmological parameters
    double h;                       // Little h = H0/(100km/s/Mpc)
    double OmegaB;                  // Baryon density today
    double OmegaCDM;                // CDM density today
    double OmegaLambda;             // Dark energy density today
    double Neff;                    // Effective number of relativistic species (3.046 or 0 if ignoring neutrinos)
    double TCMB;                    // Temperature of the CMB today in Kelvin
   
    // Derived parameters
    double OmegaR;                  // Photon density today (follows from TCMB)
    double OmegaNu;                 // Neutrino density today (follows from TCMB and Neff)
    double OmegaK;                  // Curvature density = 1 - OmegaM - OmegaR - OmegaNu - OmegaLambda
    double H0;                      // The Hubble parameter today H0 = 100h km/s/Mpc
    double OmegaM_tot;
    double OmegaRad_tot;

    // Start and end of x-integration (can be changed)
    double x_start = Constants.x_start;
    double x_end   = 5.0;//Constants.x_end;
    const int nx   = 1e4;
    const int npts = 1e4;
    const double dx = (double)(x_end-x_start)/nx;

    // Splines to be made
    Spline eta_of_x_spline{"eta"};
    Spline t_of_x_spline{"t"};
 
  public:

    // Constructors 
    BackgroundCosmology() = delete;
    BackgroundCosmology(
        double h, 
        double OmegaB, 
        double OmegaCDM, 
        double OmegaK,
        double Neff, 
        double TCMB
        );

    // Print some useful info about the class
    void info() const;

    // Do all the solving
    void solve();
    void solve_eta(double x_low=Constants.x_start, 
                   double x_high=Constants.x_end);
    void solve_t();

    // Output some results to file
    void output(const std::string filename) const;
    void output_dL(
      const std::string filename,
      const double x_low,
      const double x_high) const;


    // Get functions that we must implement
    double eta_of_x(double x) const;
    double H_of_x(double x) const;
    double Hp_of_x(double x) const;
    double dHpdx_of_x(double x) const;
    double ddHpddx_of_x(double x) const;
    double Hp_over_H0_squared(double x) const;
    double get_OmegaB(double x = 0.0) const; 
    double get_OmegaR(double x = 0.0) const;
    double get_OmegaRtot(double x = 0.0) const; 
    double get_OmegaNu(double x = 0.0) const;
    double get_OmegaCDM(double x = 0.0) const; 
    double get_OmegaLambda(double x = 0.0) const; 
    double get_OmegaK(double x = 0.0) const; 
    double get_OmegaMnu(double x = 0.0) const; // WHAT???
    double get_H0() const;
    double get_h() const;
    double get_rho_c0() const;
    double get_Neff() const;
    double get_TCMB(double x = 0.0) const;

    // Time 
    double get_t_of_x(double x = 0.0) const;

    // Distance measures
    double get_luminosity_distance_of_x(double x) const;
    double get_comoving_distance_of_x(double x) const;

    double get_acceleration_onset() const;
    double get_mr_equality(Vector x_array) const;
    double get_mL_equality() const;
    double get_k_eta_equals_unity(double k) const;
};

#endif
