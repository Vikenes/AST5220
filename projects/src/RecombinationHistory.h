#ifndef _RECOMBINATION_HISTORY_HEADER
#define _RECOMBINATION_HISTORY_HEADER
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include <iomanip>
#include <string>

using Vector = std::vector<double>;

class RecombinationHistory{
  private:

    // The cosmology we use
    BackgroundCosmology *cosmo = nullptr;
    
    // Helium fraction
    double Yp;

    // Xe for when to switch between Saha and Peebles
    double Xe_saha_limit;

    // The start and end points for recombination arrays (can be modified)
    const double x_start  = -15.0;
    const double x_end    = Constants.x_end;

    // Numbers of points of Xe,ne array (modify as you see fit)
    const int npts_rec_arrays = 1e4;
  
    // const double Xe_saha_limit = 1e-6;
    
    
    double c_;
    double sigma_T_; 
    double m_e_;
    double k_b_;
    double hbar_;
    double eps0_;
    double lambda_2s1s_;



    //===============================================================
    // [1] Computation of Xe (Saha and Peebles equation)
    //===============================================================
 
    // Compute Xe from the Saha equation
    std::pair<double,double> electron_fraction_from_saha_equation(double x) const;
    
    // Right hand side of the dXedx Peebles equation
    int rhs_peebles_ode(double x, const double *y, double *dydx);
    
    // Solve for Xe 
    void solve_number_density_electrons();
    
    //===============================================================
    // [2] Compute tau and visibility functions
    //===============================================================

    // Number of points of tau, g_tilde arrays, and their derivatives
    const int npts_tau = 1e5;

    // The two things we need to solve: Xe/ne and tau
    void solve_for_optical_depth_tau();


    //===============================================================
    // [3] Compute sound horizon
    //===============================================================
    void solve_sound_horizon();
    
    
    // Splines contained in this class
    Spline Xe_of_x_spline{"Xe"};
    Spline log_ne_of_x_spline{"ne"};
    Spline tau_of_x_spline{"tau"}; 
    Spline dtau_dx_spline{"dtau_dx"};
    Spline ddtau_ddx_spline{"ddtau_dxx"};
    Spline g_tilde_of_x_spline{"g"};  
    Spline dg_dx_spline{"dg_dx"};
    Spline ddg_ddx_spline{"ddg_dxx"};
    Spline s_of_x_spline{"s(x)"};

    // Number of points written to file
    const int nx_write = 5e4;

  public:

    // Construtors
    RecombinationHistory() = delete;
    RecombinationHistory(
        BackgroundCosmology *cosmo, 
        double Yp,
        double Xe_saha_limit);

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Output some data to file
    void output(const std::string filename) const;
    void output_important_times(const std::string filename) const;


    // Get functions that we must implement
    double tau_of_x(double x) const;
    double dtaudx_of_x(double x) const;
    double ddtauddx_of_x(double x) const;
    double g_tilde_of_x(double x) const;
    double dgdx_tilde_of_x(double x) const;
    double ddgddx_tilde_of_x(double x) const;
    double s_of_x(double x) const;
    double nb_of_x(double x) const;
    double Xe_of_x(double x) const;
    double ne_of_x(double x) const;
    double cs_of_x(double x) const;
    double get_Yp() const;

    double get_x_at_decoupling_tau(
        double xmin=Constants.x_start,
        double xmax=Constants.x_end) const;
        
    double get_x_at_decoupling_g(
        double xmin=Constants.x_start,
        double xmax=Constants.x_end) const;
    double get_x_at_recombination(double Xe_rec=0.1) const;

};

#endif
