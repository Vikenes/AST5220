#ifndef _PERTURBATIONS_HEADER
#define _PERTURBATIONS_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <vector>
#include <fstream>
#include <algorithm>
#include <filesystem>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class Perturbations{
  private:

    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;
   
    // The scales we integrate over
    const int n_k        = 200;
    const double k_min   = Constants.k_min;
    const double k_max   = Constants.k_max;
    
    // Start and end of the time-integration
    const int n_x        = 5000;
    const double x_start = Constants.x_start;
    const double x_end   = Constants.x_end;

    Vector x_array_full = Utils::linspace(x_start, x_end, n_x);

    // Below is a full list of splines you probably need, 
    // but you only need to make the splines you will need

    // Splines of scalar perturbations quantities
    Spline2D delta_cdm_spline{"delta_cdm_spline"};
    Spline2D delta_b_spline{"delta_b_spline"};
    Spline2D v_cdm_spline{"v_cdm_spline"};
    Spline2D v_b_spline{"v_b_spline"};
    Spline2D Phi_spline{"Phi_spline"};
    Spline2D Pi_spline{"Pi_spline"};
    Spline2D Psi_spline{"Psi_spline"};
   
    // Splines of source functions (ST for temperature; SE for polarization)
    Spline2D ST_spline{"ST_spline"};
    std::vector<Spline2D> ST_component_splines;
    
    // Splines of mulipole quantities
    std::vector<Spline2D> Theta_splines;

    double c_;
    double H0_; 
    double H0_squared_;
    double OmegaR0_;
    double OmegaB0_;
    double OmegaCDM0_;


    //==========================================================
    // [1] Tight coupling ODE system
    //==========================================================

    // Set the initial conditions at the start (which is in tight coupling)
    Vector set_ic(
        const double x, 
        const double k) const;
    
    // Right hand side of the ODE in the tight coupling regime
    int rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx);

    // Compute the time when tight coupling ends
    std::pair<double,int> get_tight_coupling_time(const double k) const;

    //==========================================================
    // [2] The full ODE system 
    //==========================================================
  
    // Set initial condition after tight coupling
    Vector set_ic_after_tight_coupling(
        const Vector &y_tight_coupling, 
        const double x, 
        const double k) const;

    // Right hand side of the ODE in the full regime
    int rhs_full_ode(double x, double k, const double *y, double *dydx);
    
    //==========================================================
    // [3] Integrate the full system
    //==========================================================
    void integrate_perturbations();
    
    //==========================================================
    // [4] Compute source functions from the result
    //==========================================================
    void compute_source_functions(int term=0);

  public:

    // Constructors
    Perturbations() = default;
    Perturbations(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec); 

    // Do all the solving
    void solve(bool source=true, int term=0);
    
    // Print some useful info about the class
    void info() const;

    // Output info to file
    bool check_existence(std::string filename);
    void output(const double k, const std::string filename) const;
    void outputTheta0(const Vector &k_arr, std::string filename);
    void outputPsi(const Vector &k_arr, std::string filename);


    // Get the quantities we have integrated
    double get_delta_cdm(const double x, const double k) const;
    double get_delta_b(const double x, const double k) const;
    double get_v_cdm(const double x, const double k) const;
    double get_v_b(const double x, const double k) const;
    double get_Phi(const double x, const double k) const;
    double get_Psi(const double x, const double k) const;
    double get_Pi(const double x, const double k) const;
    double get_Theta(const double x, const double k, const int ell) const;

    double compute_Theta2_tc(const double x, const double k, const double Theta1) const;
    double compute_Psi(const double x, const double k, const double Theta2, const double Phi) const;
    double get_Source_T(const double x, const double k) const;
    double get_Source_T_component(const double x, const double k, const int term) const;

};

#endif
