#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{
  
}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  const double eta0 = cosmo->eta_of_x(0.0);
  const double dk = 2.0 * M_PI / eta0 / (double)los_samples_per_osc;
  const int nk = int((k_max - k_min) / dk);

  const double dlogk = 2.0 * M_PI / eta0 / (double)cell_samples_per_osc;
  const int nlogk = int((k_max - k_min) / dlogk);

  Vector k_array = Utils::linspace(k_min, k_max, nk);
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), nlogk);

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  

}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
  
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================

  const double zmax = k_max * cosmo->eta_of_x(0.0);
  const double zmin = 0.0; 
  // const int n_samples_per_oscillation = 20;
  const double dz = 2.0 * M_PI / (double)bessel_samples_per_osc; 
  const int n_z = (zmax - zmin) / dz;

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];

    Vector z_array = Utils::linspace(zmin, zmax, n_z);
    Vector j_ell_array(n_z);
    for(int iz=0; iz < n_z; iz++){
      j_ell_array[iz] = Utils::j_ell(ell, z_array[iz]);
    }
    // std::string j_ell_name = "j_" + std::to_string(ell); 
    j_ell_splines[i].create(z_array, j_ell_array);//, j_ell_name);
  }
  

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  const int nells = ells.size();
  const double eta0 = cosmo->eta_of_x(0.0);
  const double dx = (x_end - x_start) / (double)n_x;

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t ik = 0; ik < k_array.size(); ik++){

    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    double k = k_array[ik];
    for(int iell=0; iell<nells; iell++){
      std::function<double(double)> integrand = [&](double x){
        return source_function(x, k) * j_ell_splines[iell](k*(eta0 - cosmo->eta_of_x(x)));
      };
      result[iell][ik] = integrate_trapezoidal(integrand, x_start, x_end, dx); 
    }

    // Store the result for Source_ell(k) in results[ell][ik]
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);
  // Spline the result and store it in thetaT_ell_of_k_spline
  // ...
  for(int i=0; i<nells; i++){
    int ell = ells[i];
    // std::string Theta_ell_name = "Theta_" + std::to_string(ell) + "of_k_spline";
    thetaT_ell_of_k_spline[i].create(k_array, thetaT_ell_of_k[i]);//, Theta_ell_name);
  }

}



//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();
  
  const double log_k_min = log_k_array[0];
  const double log_k_max = log_k_array[log_k_array.size()-1];
  const double dlogk = log_k_array[1] - log_k_array[0];

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================
  Vector result(nells);
  Utils::StartTiming("C_ell");
  for(int iell=0; iell<nells; iell++){

    std::function<double(double)> integrand = [&](double logk){
          double k = exp(logk);
          return primordial_power_spectrum(k) 
                  * f_ell_spline[iell](k) 
                  * g_ell_spline[iell](k);
        };
    result[iell] = 4.0 * M_PI * integrate_trapezoidal(integrand, 
                                                      log_k_min, 
                                                      log_k_max, 
                                                      dlogk);
  }
  Utils::EndTiming("C_ell");


  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================

  // ...
  // ...
  // ...

  return pofk;
}


double PowerSpectrum::integrate_trapezoidal(
    std::function<double(double)> &f, 
    const double zmin,
    const double zmax,
    const double dz){

  double sum = 0.0;
  double z = zmin;

  sum += 0.5 * f(z);
  z += dz;

  while(z < zmax){
    sum += f(z);
    z += dz;
  }
  z = zmax;
  sum += 0.5 * f(z);

  return sum * dz;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };

  std::cout << "Writing file to: " << filename << std::endl;

  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

