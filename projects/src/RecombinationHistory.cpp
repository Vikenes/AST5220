#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();

  // Compute and spline s(x)
  solve_sound_horizon();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================
void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  

  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // Get X_e from solving the Saha equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit){
      saha_regime = false;
    }

      

    if(saha_regime){
      //=============================================================================
      // Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;


    } else {
      //==============================================================
      // Compute X_e from current time till today by solving 
      // the Peebles equation
      //==============================================================

      std::cout << "End of Saha regime reached at x=" << x_array[i-1] << std::endl; 

      // The Peebles ODE equation
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };

      // Final value of i where Saha holds  
      int saha_end_ix = i-1;

      // Create x-array for the peebles regime   
      Vector x_array_peebles = Vector(x_array.begin()+saha_end_ix, x_array.end()); 
      
      // Use final Saha value as IC for peebles
      Vector Xe_ic{Xe_arr[saha_end_ix]};

      // Solve ODE
      ODESolver peebles_Xe_ode;
      peebles_Xe_ode.solve(dXedx, x_array_peebles, Xe_ic);

      auto Xe_arr_peebles = peebles_Xe_ode.get_data_by_component(0);

      // Store the result we get from the Peebles equation  
      for(int i=saha_end_ix; i < npts_rec_arrays; i++){
        Xe_arr[i] = Xe_arr_peebles[i - saha_end_ix];
        ne_arr[i] = Xe_arr_peebles[i - saha_end_ix] * nb_of_x(x_array[i]); 
      }

      // Computation finished after completing Peebles 
      break;
    }
  }

  //=============================================================================
  // Splining the result.  
  //=============================================================================
  Vector log_ne_arr = log(ne_arr);

  Xe_of_x_spline.create(x_array, Xe_arr, "log_Xe");
  log_ne_of_x_spline.create(x_array, log_ne_arr, "log_ne");

  Utils::EndTiming("Xe");
  
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  // Physical constants
  double k_b         = Constants.k_b;
  double m_e         = Constants.m_e;
  double hbar        = Constants.hbar;
  double epsilon_0   = Constants.epsilon_0;

  // Fetch cosmological parameters
  const double Tb          = cosmo->get_TCMB(x);
  const double kT          = k_b * Tb;
  const double nb          = nb_of_x(x); 


  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================

  double pre_factor      = m_e * k_b*Tb / (2.0 * M_PI * hbar * hbar);
  double y               = 1.0/nb * pow(pre_factor, 1.5) * exp(-epsilon_0 / (kT));

  // Compute Xe 
  // Set Xe=1 when y is very large to avoid numerical errors. 
  if(y > 1e7){
    Xe = 1.0;
  } 
  else {
    Xe =  0.5*y * (sqrt(1.0 + 4.0 / y) - 1.0);  
  };

  // Compute electron density 
  ne = Xe * nb; 

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of X_e
  double X_e         = Xe[0];

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  const double Tb           = cosmo->get_TCMB(x);
  const double H            = cosmo->H_of_x(x);
  const double nb           = nb_of_x(x); // Assume all baryons are protons 

  // Relevant physical quantities for the Peebles equation
  double kT                 = k_b * Tb;
  double eps0_over_kT       = epsilon_0 / kT;
  double eps_over_hc_cubed  = epsilon_0*epsilon_0*epsilon_0 / (hbar*hbar*hbar * c*c*c);

  //=============================================================================
  // Compute parameters of the Peebles equation 
  //=============================================================================
  double phi2      = 0.448 * log(eps0_over_kT);
  double alpha2    = 8. / sqrt(3.*M_PI) * c * sigma_T * sqrt(eps0_over_kT) * phi2;
  double beta      = alpha2 * pow((m_e*kT / (2.*M_PI*hbar*hbar)),1.5) * exp(-eps0_over_kT);
  double nH        = (1. - Yp) * nb; // Yp=0. 
  double n1s       = (1. - X_e) * nH;

  // Compute beta2. Avoid overflow at low temperature 
  double beta2 = 0.0;
  if(eps0_over_kT<200)
    beta2 = beta * exp(0.75 * eps0_over_kT);

  double Lambda_alpha = H * 27. * eps_over_hc_cubed / (64. * M_PI*M_PI * n1s);
  double Cr           = (lambda_2s1s + Lambda_alpha) 
                         / (lambda_2s1s + Lambda_alpha + beta2); 

  dXedx[0] = Cr/H * (beta*(1. - X_e) - nH*alpha2 * X_e*X_e);

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================
void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  Vector x_array = Utils::linspace(x_start, 0, npts_tau);
  Vector x_backwards_array = Utils::linspace(0, x_start, npts_tau);


  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    double c       = Constants.c; 
    double sigma_T = Constants.sigma_T;
    dtaudx[0] = - c * sigma_T * ne_of_x(x) * exp(x) / cosmo->Hp_of_x(x); 
    return GSL_SUCCESS;
  };

  // Solve backwards, starting from x=0
  Vector tau_ic_back{0.0}; // tau(x=0) = 0   

  // Solve ODE 
  ODESolver tau_ode_backwards;
  tau_ode_backwards.set_accuracy(1e-5, 1e-8, 1e-8);

  tau_ode_backwards.solve(dtaudx, x_backwards_array, tau_ic_back);

  // 
  auto tau_array_backwards      = tau_ode_backwards.get_data_by_component(0);
  auto dtau_dx_array_backwards  = tau_ode_backwards.get_derivative_data_by_component(0);
  
  // Store tau for increasing x 
  // Compute analytical expressions for t'(x) and g_tilde(x) 
  Vector tau_array(npts_tau);
  Vector dtau_dx_array(npts_tau);

  std::reverse(tau_array_backwards.begin(), tau_array_backwards.end());
  std::reverse(dtau_dx_array_backwards.begin(), dtau_dx_array_backwards.end());

  std::copy(tau_array_backwards.begin(), tau_array_backwards.end(), tau_array.begin());
  std::copy(dtau_dx_array_backwards.begin(), dtau_dx_array_backwards.end(), dtau_dx_array.begin());

  // Vector dtau_dx_array(npts_tau);
  Vector g_tilde_array(npts_tau);

  for(int i=0; i<npts_tau; i++){
    double x_ = x_array[i];

    // tau_array[i] = tau_array_backwards[npts_tau-i-1];
    // dtau_dx_array[i] = dtau_dx_array_backwards[npts_tau-i-1];
    // dtau_dx_array[i] = - ne_of_x(x_) * Constants.c * Constants.sigma_T / cosmo->H_of_x(x_);
    g_tilde_array[i] = - dtau_dx_array[i] * exp(-tau_array[i]);
  }

  // Spline results 
  tau_of_x_spline.create(x_array, tau_array, "tau_of_x");
  dtau_dx_spline.create(x_array, dtau_dx_array, "dtau_dx");
  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g_tilde_of_x");

  Vector ddtau_ddx_array(npts_tau);
  for (int i=0; i<npts_tau; i++){
    ddtau_ddx_array[i] = dtau_dx_spline.deriv_x(x_array[i]);
  }  

  ddtau_ddx_spline.create(x_array, ddtau_ddx_array, "ddtau");

  // Create spline for g'(x). 
  // Avoid numerical errors when differentiating g twice.   
  Vector dg_dx_array(npts_tau);
  
  for(int i=0; i<npts_tau; i++){
    double x_ = x_array[i];
    dg_dx_array[i] = (-ddtauddx_of_x(x_) + dtaudx_of_x(x_)*dtaudx_of_x(x_))* exp(-tau_of_x(x_));
  }

  dg_dx_spline.create(x_array, dg_dx_array, "dg_dx");

  Utils::EndTiming("opticaldepth");
}


//====================================================
// Solve the sound horizon from the cosmology. Spline result
//====================================================
void RecombinationHistory::solve_sound_horizon(){
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);

  // The ODE for ds/dx
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){
    dsdx[0] = cs_of_x(x) / cosmo->Hp_of_x(x);
    return GSL_SUCCESS;
  };
  
  // The initial condition for ds/dx
  double s_init = cs_of_x(x_start) / cosmo->Hp_of_x(x_start);
  Vector s_ic{s_init};

  // Solve the ODE
  ODESolver s_ode;
  s_ode.solve(dsdx, x_array, s_ic);

  // Spline result 
  auto s_array = s_ode.get_data_by_component(0);
  s_of_x_spline.create(x_array, s_array, "s(x)"); 

}


//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return dtau_dx_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  // return dtau_dx_spline.deriv_x(x);
  return ddtau_ddx_spline(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return dg_dx_spline(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return dg_dx_spline.deriv_x(x);
}

double RecombinationHistory::s_of_x(double x) const{
  return s_of_x_spline(x);
}

double RecombinationHistory::nb_of_x(double x) const{
  //=============================================================================
  // Compute baryon number density at a given x
  //=============================================================================
  double a_cubed_inv = exp(-3.0*x);
  double OmegaB      = cosmo->get_OmegaB(); 
  double rho_c0      = cosmo->get_rho_c0();
  double nb          = OmegaB * rho_c0  / Constants.m_H / exp(3.*x);  
  return nb;
}


double RecombinationHistory::Xe_of_x(double x) const{
  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::cs_of_x(double x) const{
  //====================================================
  // Compute sound speed of photon baryon plasma
  //====================================================
  double c            = Constants.c;
  double OmegaR       = cosmo->get_OmegaR(x);
  double OmegaB       = cosmo->get_OmegaB(x);
  double nominator    = 4.0*OmegaR;
  double denominator  = 3.0*OmegaB + nominator;
  double c_s          = c / sqrt(3.0) * sqrt(nominator / denominator);

  return c_s;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  
  std::string fname = filename;
  if(Xe_saha_limit < 0.99){
    int idx = filename.find(".");
    fname = fname.insert(idx, "_saha");
  }

  std::ofstream fp(fname.c_str());

  std::cout << "Writing result to " << fname << std::endl; 

  const int npts       = nx_write;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  fp << "nx=" << npts << ", npts_rec=" << npts_rec_arrays << ", npts_tau_g=" << npts_tau << "\n";
  fp << "x Xe ne tau dtaudx ddtauddx g dgdx ddgddx s \n";
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << s_of_x(x)            << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

void RecombinationHistory::output_important_times(const std::string filename) const{
  std::string fname = filename;

  if(Xe_saha_limit < 0.99){
    int idx = filename.find(".");
    fname = fname.insert(idx, "_saha");
  }

  std::cout << fname << std::endl;
  
  std::ofstream fp(fname.c_str());

  std::cout << "Computing decoupling and recombination time. " << std::endl
  << "Writing to: " << fname << std::endl;

  // Locate tau=1. 
  double tau_at_dec = 1.0; 
  auto x_range_tau = std::pair<double,double>(x_start, x_end);
  double x_dec_tau = Utils::binary_search_for_value(tau_of_x_spline, 
                                                  tau_at_dec,
                                                  x_range_tau);

  // Find peak of visibility function 
  double dg_dx_at_dec = 0.0;
  auto x_range_g_peak = std::pair<double,double>(x_dec_tau - 0.1, x_dec_tau + 0.1);
  double x_dec_g_peak = Utils::binary_search_for_value(dg_dx_spline, 
                                                      dg_dx_at_dec,
                                                      x_range_g_peak);
  
  // Find recombination time 
  double log_Xe_at_recomb = 0.1;
  auto x_range_recomb = std::pair<double,double>(x_start, x_end);
  double x_recomb = Utils::binary_search_for_value(Xe_of_x_spline,
                                                      log_Xe_at_recomb,
                                                      x_range_recomb);

  
  cosmo->solve_t();
  auto z = [&](double x){ return exp(-x) - 1.0; };
  auto t = [&](double x){ return cosmo->get_t_of_x(x); }; 

  fp << "LS Recombination \n";
  fp << "x: " << x_dec_g_peak << " " << x_recomb << " \n";
  fp << "z: " << z(x_dec_g_peak) << " " << z(x_recomb) << " \n"; 
  fp << "t: " << t(x_dec_g_peak) << " " << t(x_recomb) << " \n"; 
  fp << "r: " << s_of_x(x_dec_g_peak) << " " << s_of_x(x_recomb) << " \n"; 



  // Spline GH_saha{"gh"};
  // Spline GH_peeb{"gh"};
  // Vector x = Utils::linspace(-10, -5, npts_rec_arrays);
  // Vector GGH_saha(npts_rec_arrays);
  // Vector GGH_peeb(npts_rec_arrays);
  // double sigmaT = Constants.sigma_T;
  // double c = Constants.c;
  // for(int i=0; i<npts_rec_arrays; i++){
  //   auto Xe_ne_saha = electron_fraction_from_saha_equation(x[i]);
  //   double ne_saha = Xe_ne_saha.second;
  //   double ne_peeb = ne_of_x(x[i]);
  //   double H = cosmo->H_of_x(x[i]);
  //   GGH_saha[i] = sigmaT*c * ne_saha - H;
  //   GGH_peeb[i] = sigmaT*c * ne_peeb - H;
  // }
  // GH_saha.create(x, GGH_saha, "gghh_saha");
  // GH_peeb.create(x, GGH_peeb, "gghh_peeb");


  // double xdec_saha = Utils::binary_search_for_value(GH_saha, 0.0, std::pair<double,double>(-10, -5));
  // double xdec_peeb = Utils::binary_search_for_value(GH_peeb, 0.0, std::pair<double,double>(-10, -5));

  // std::cout << "gamma=H at" << std::endl;
  // std::cout << "Saha: x=" << std::setprecision(10) << xdec_saha << ", z=" << std::setprecision(10) << z_of_x(xdec_saha) << std::endl;
  // std::cout << "Peeb: x=" << std::setprecision(10) << xdec_peeb << ", z=" << std::setprecision(10) << z_of_x(xdec_peeb) << std::endl;


  /*
  double z_dec_tau       = z_of_x(x_dec_tau);
  double z_dec_g_peak    = z_of_x(x_dec_g_peak);
  double z_recomb        = z_of_x(x_recomb);

  double t_dec_tau       = cosmo->get_t_of_x(x_dec_tau);
  double t_dec_g_peak    = cosmo->get_t_of_x(x_dec_tau);
  double t_recomb        = cosmo->get_t_of_x(x_recomb);

  double r_s_dec_tau     = s_of_x(x_dec_tau);
  double r_s_dec_g_peak  = s_of_x(x_dec_g_peak);
  double r_s_recomb      = s_of_x(x_recomb);



  std::cout << "Decoupling: " << "tau = " << tau_of_x(x_dec_tau) << std::endl;
  std::cout << "   x  = " << x_dec_tau << std::endl;
  std::cout << "   z  = " << z_dec_tau << std::endl;
  std::cout << "   t  = " << t_dec_tau << std::endl;
  std::cout << "   rs = " << r_s_dec_tau << std::endl << std::endl;

  std::cout << "Decoupling: " << "gmax = " << g_tilde_of_x(x_dec_g_peak)
  << ", dgdx = " << dgdx_tilde_of_x(x_dec_g_peak) << std::endl;
  std::cout << "   x  = " << x_dec_g_peak << std::endl;
  std::cout << "   z  = " << z_dec_g_peak << std::endl;
  std::cout << "   t  = " << t_dec_g_peak << std::endl;
  std::cout << "   rs = " << r_s_dec_g_peak << std::endl << std::endl;

  std::cout << "Recombination: " << "Xe = " << Xe_of_x(x_recomb) << std::endl;
  std::cout << "   x  = " << x_recomb << std::endl;
  std::cout << "   z  = " << z_recomb << std::endl;
  std::cout << "   t  = " << t_recomb << std::endl;
  std::cout << "   rs = " << r_s_recomb << std::endl;

  */



}