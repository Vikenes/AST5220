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
  c_                = Constants.c;
  eta0_             = cosmo->eta_of_x(0.0);
  TWO_PI_over_eta0_ = 2.0 * M_PI / eta0_;
  nells_            = ells.size();
  H0_               = cosmo->get_H0();
  OmegaMtot_        = cosmo->get_OmegaB() + cosmo->get_OmegaCDM();
  Mpc_              = Constants.Mpc;
}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(bool load_data){
  
  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  // const double eta0 = cosmo->eta_of_x(0.0);

  // const double dk      = k_stepsize_from_N_osc_samples(los_samples_per_osc);
  const int nk         = n_k_from_N_osc_samples(los_samples_per_osc);
  // const double dlogk   = k_stepsize_from_N_osc_samples(cell_samples_per_osc);
  const int nlogk      = n_k_from_N_osc_samples(cell_samples_per_osc);

  Vector k_array = Utils::linspace(k_min, k_max, nk);
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), nlogk);

  
  //=========================================================================
  // Load data from previous run.
  // Mainly for testing purposes to avoid long simulations. 
  //=========================================================================


  //=========================================================================
  // Generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // Line of sight integration to get Theta_ell(k)
  //=========================================================================
  line_of_sight_integration(k_array);

  //=========================================================================
  // Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
}

double PowerSpectrum::k_stepsize_from_N_osc_samples(int samples_per_osc) const{
  return TWO_PI_over_eta0_ / (double)samples_per_osc;
}

double PowerSpectrum::n_k_from_N_osc_samples(int samples_per_osc) const{
  return int((k_max - k_min) / k_stepsize_from_N_osc_samples(samples_per_osc));
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

  const double zmax = k_max * eta0_;
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
  // const int nells = ells.size();
  const double dx = (x_end - x_start) / (double)n_x;

  Vector x_LOS_array = Utils::linspace(x_start_LOS, x_end, n_x);

  // const int n_x_LOS = int((x_end - x_start_LOS) / dx);
  // Vector x_LOS = Utils::linspace(x_start_LOS, x_end, n_x_LOS);

  // std::cout << "nxlos=" << n_x_LOS << std::endl;
  // return result;   
  // Vector x_LOS;
  // for(int ix=0; ix<n_x; ix++){
  //   double x_ = x_array[ix];
  //   if(rec->g_tilde_of_x(x_) > 1e-5 || x_ > -4.){
  //     x_LOS.push_back(x_);
  //   }
  // }
  // std::cout << "N=" << x_LOS.size() << std::endl;
  // std::cout << "x_0=" << x_LOS[0] << std::endl;


  for(size_t ik = 0; ik < k_array.size(); ik++){

    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    double k = k_array[ik];
   
    #pragma omp parallel for schedule(dynamic, 1)
    for(int iell=0; iell<nells_; iell++){
      // std::function<double(double)> integrand = [&](double x){
        // return source_function(x, k) * j_ell_splines[iell](k*(eta0 - cosmo->eta_of_x(x)));
      // };
      double integral_sum = 0;
      // double x;  
      double x = x_start_LOS;  
      for(int ix=0; ix<x_LOS_array.size(); ix++){
        // x = x_LOS[ix];
        integral_sum += source_function(x, k) * j_ell_splines[iell](k*(eta0_ - cosmo->eta_of_x(x)));
        x += dx;
      }
      result[iell][ik] = integral_sum * dx; 
    }

    // Store the result for Source_ell(k) in results[ell][ik]
  }
  std::cout << std::endl;

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  // const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells_);

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
  for(int i=0; i<nells_; i++){
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
  // const int nells      = ells.size();
  
  const double log_k_min = log_k_array[0];
  const double log_k_max = log_k_array[log_k_array.size()-1];
  const double dlogk = log_k_array[1] - log_k_array[0];

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================
  Vector result(nells_);
  Utils::StartTiming("C_ell");
  for(int iell=0; iell<nells_; iell++){

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
  return A_s * pow( Mpc_ * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================
  double k_SI = k_mpc / Mpc_;
  double Phi = pert->get_Phi(x, k_SI);
  double a = 1.;

  double Delta_M = c_*c_*k_SI*k_SI * Phi / (3./2. * OmegaMtot_ * H0_*H0_);
  double k_mpc_cubed = k_mpc * k_mpc * k_mpc;
  double P_primordial = 2.0*M_PI*M_PI * primordial_power_spectrum(k_SI) / k_mpc_cubed; 

  pofk = abs(Delta_M*Delta_M) * P_primordial;
  return pofk;
}

double PowerSpectrum::get_k_eq() const{
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  double x_eq = cosmo->get_mr_equality(x_array);
  double a_eq = exp(x_eq);
  double k_eq = a_eq * cosmo->H_of_x(x_eq) / c_;

  return k_eq;
}

double PowerSpectrum::get_Theta_squared_over_k(const int iell, const double k) const{
  double Theta = thetaT_ell_of_k_spline[iell](k);
  return abs(Theta*Theta) / k; 
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
  int nk    = n_k_from_N_osc_samples(los_samples_per_osc);
  int nlogk = n_k_from_N_osc_samples(cell_samples_per_osc);
  int pos = filename.find(".txt");
  if (pos != std::string::npos) {
    std::string file_info;
    file_info  = "_nx" + std::to_string(n_x);
    file_info += "_nk" + std::to_string(nk);
    file_info += "_nlogk" + std::to_string(nlogk);
    filename.insert(pos, file_info);
  }

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


void PowerSpectrum::outputThetas(std::string filename, int nk) const{

  int pos = filename.find(".txt");

  if (pos != std::string::npos) {
    std::string file_info  = "_nk" + std::to_string(nk) + "_nx" + std::to_string(n_x);
    filename.insert(pos, file_info);
  }

  std::ofstream fp(filename.c_str());


  std::string header = "k*eta0, ell= ";
  for(int iell=0; iell<nells_; iell++){
    header += std::to_string(int(ells[iell])) + " ";
  }
  fp << "eta0= " << eta0_ << "\n";
  fp << header << "\n";

  auto k_log_values = Utils::linspace(log(k_min), log(k_max), nk);
  auto k_values = exp(k_log_values);

  auto print_data = [&] (const double k) {
    fp << k                           << " ";
    for(int iell=0; iell<nells_; iell++){
      fp << thetaT_ell_of_k_spline[iell](k) << " ";
    }
    fp << "\n";
  };

  std::cout << "Writing file to: " << filename << std::endl;
  std::for_each(k_values.begin(), k_values.end(), print_data);
  std::cout << "Finished writing to file" << std::endl;

}

void PowerSpectrum::outputCellIntegrand(std::string filename, int nk) const{

  int pos = filename.find(".txt");

  if (pos != std::string::npos) {
    std::string file_info  = "_nk" + std::to_string(nk) + "_nx" + std::to_string(n_x);
    filename.insert(pos, file_info);
  }

  std::ofstream fp(filename.c_str());

  auto kvalues = Utils::linspace(k_min, k_max, nk);

  std::string header = "k,ell: ";
  for(int iell=0; iell<nells_; iell++){
    header += std::to_string(int(ells[iell])) + " ";
  }
  fp << header << "\n";

  auto print_data = [&] (const double k) {
    fp << k                                     << " ";
    for(int iell=0; iell<nells_; iell++){
      fp << get_Theta_squared_over_k(iell, k) << " ";
    }
    fp << "\n";
  };

  std::cout << "Writing file to: " << filename << std::endl;
  std::for_each(kvalues.begin(), kvalues.end(), print_data);
  std::cout << "Finished writing to file" << std::endl;

}


void PowerSpectrum::outputPS(std::string filename, int nk) const{
  // Output in standard units of muK^2
  int pos = filename.find(".txt");
  if (pos != std::string::npos) {
    std::string file_info  = "_nk" + std::to_string(nk);
    filename.insert(pos, file_info);
  }

  std::ofstream fp(filename.c_str());
  
  auto kvalues = Utils::linspace(k_min*Mpc_, k_max*Mpc_, nk);

  const double h = cosmo->get_h();
  const double knorm = 1.0 / h;
  const double PSnorm = h*h*h; 


  auto print_data = [&] (const double k) {
    double k_h_Mpc = k * knorm;
    fp << k_h_Mpc                                     << " ";
    fp << get_matter_power_spectrum(0.0, k) * PSnorm  << " ";
    fp << "\n";
  };

  double k_eq = get_k_eq();

  std::cout << "Writing file to: " << filename << std::endl;
  fp << k_eq * Mpc_ / h << "\n";
  std::for_each(kvalues.begin(), kvalues.end(), print_data);

  std::cout << "Finished writing to file" << std::endl;

}