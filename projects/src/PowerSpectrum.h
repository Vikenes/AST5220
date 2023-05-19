#ifndef _POWERSPECTRUM_HEADER
#define _POWERSPECTRUM_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <functional>
#include <utility> 
#include <fstream> 
#include <algorithm>
#include <filesystem>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class PowerSpectrum {
  private:

    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;
    Perturbations *pert        = nullptr;

    // Parameters defining the primordial power-spectrum
    double A_s        = 2.1e-9;
    double n_s        = 0.965;
    double kpivot_mpc = 0.05;

    // The k-values we compute Theta_ell(k) etc. for
    // const int nk      = 100;
    const double k_min = Constants.k_min;
    const double k_max = Constants.k_max;

    const int los_samples_per_osc = 16;
    const int bessel_samples_per_osc = 20;
    const int cell_samples_per_osc = 16;

    // The x-values we integrate to compute Theta_ell(k) etc. for
    const int n_x        = 1500;
    const double x_start = -15.0;
    const double x_end   = Constants.x_end;
    const double x_start_LOS = x_start;// = -11.0;

    
    // The ells's we will compute Theta_ell and Cell for
    Vector ells{ 
        2,    3,    4,    5,    6,    7,    8,    10,   12,   15,   
        20,   25,   30,   40,   50,   60,   70,   80,   90,   100,  
        120,  140,  160,  180,  200,  225,  250,  275,  300,  350,  
        400,  450,  500,  550,  600,  650,  700,  750,  800,  850,  
        900,  950,  1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 
        1900, 1950, 2000};
   
    // For initial testing 
    // Vector ells{ 
    //     2,    5,    10,     20,   
    //     35,   50,   100,    160,  
    //     225,  300,  400,    550, 
    //     700,  850,  1000,   1200, 
    //     1400, 1550, 1700,   1850, 
    //     2000};


    // Constants being used throughout 
    double c_;
    double eta0_;
    int nells_;
    double H0_;
    double OmegaMtot_;
    double TWO_PI_over_eta0_;
    double Mpc_;

    //=====================================================================
    // [1] Create bessel function splines needed for the LOS integration
    //=====================================================================

    // Splines of bessel-functions for each value of ell in the array above
    std::vector<Spline> j_ell_splines;
    
    // Generate splines of bessel-functions for each ell needed
    // to do the LOS integration
    void generate_bessel_function_splines();
    
    //=====================================================================
    // [2] Do the line of sight integration and spline the result
    //=====================================================================
    
    // Do LOS integration for all ells and all k's in the given k_array
    // and for all the source functions (temperature, polarization, ...)
    void line_of_sight_integration(Vector & k_array);
  
    // Do the line of sight integration for a single quantity
    // for all ells by providing a source_function(x,k) (can be temp, pol, ...)
    Vector2D line_of_sight_integration_single(
        Vector & k_array, 
        std::function<double(double,double)> &source_function);
    
    // Splines of the reusult of the LOS integration
    std::vector<Spline> thetaT_ell_of_k_spline;
    std::vector<Spline> thetaT_ell_of_k_comp_spline;
    
    
    //=====================================================================
    // [3] Integrate to get power-spectrum
    //=====================================================================
    
    // General method to solve for Cells (allowing for cross-correlations)
    // For auto spectrum (C_TT) then call with f_ell = g_ell = theta_ell
    // For polarization C_TE call with f_ell = theta_ell and g_ell = thetaE_ell
    Vector solve_for_cell(
        Vector & logk_array,
        std::vector<Spline> & f_ell, 
        std::vector<Spline> & g_ell);

    // Splines with the power-spectra
    Spline cell_TT_spline{"cell_TT_spline"};
    Spline cell_TE_spline{"cell_TE_spline"};
    Spline cell_EE_spline{"cell_EE_spline"};

    std::vector<Spline> cell_TT_component_spline;

  public:

    // Constructors
    PowerSpectrum() = delete;
    PowerSpectrum(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec, 
        Perturbations *pert,
        double A_s,
        double n_s,
        double kpivot_mpc);
    
    // Do all the solving: bessel functions, LOS integration and then compute Cells
    void solve();
    void solve_components();

    // The dimensionless primordial power-spectrum Delta = 2pi^2/k^3 P(k)
    double primordial_power_spectrum(const double k) const;

    // Get P(k,x) for a given x in units of (Mpc)^3
    double get_matter_power_spectrum(const double x, const double k_mpc) const;

    double get_k_eq() const;
    double get_Theta_squared_over_k(const int iell, const double k) const;

    double k_stepsize_from_N_osc_samples(int samples_per_osc) const;
    double n_k_from_N_osc_samples(int samples_per_osc) const;
    double integrate_trapezoidal(
        std::function<double(double)> &f, 
        const double zmin,
        const double zmax,
        const double dz);

    // Get the quantities we have computed
    double get_cell_TT(const double ell) const;
    double get_cell_TT_component(const double ell, const int term) const;
    double get_cell_TE(const double ell) const;
    double get_cell_EE(const double ell) const;

    void add_file_info(std::string &filename, int nx);
    void add_file_info(std::string &filename);
    bool check_existence(std::string filename);

    // Output Cells in units of l(l+1)/2pi (muK)^2
    void output(std::string filename);
    void output_Cell_components(std::string filename);
    void outputThetas(std::string filename, int nk);
    void outputPS(std::string filename, int nk);

};

#endif
