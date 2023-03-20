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
    if(Xe_current <= Xe_saha_limit){
      std::cout << "Xe_current=" << Xe_current << ", x=" << x_array[i] << std::endl;
      // for(int nn=0; nn<i; nn++){
        // std::cout << "x=" << x_array[nn] << ", Xe=" << Xe_arr[nn] << std::endl;
      // }
      saha_regime = false;
    }

      

    if(saha_regime){
      //=============================================================================
      // Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;


    } else {
      std::cout << "Peebles. Xe=" << Xe_current << "/" << Xe_arr[i-1] 
      << ", xi-1=" << x_array[i-1] << ", xi=" << x_array[i] << std::endl; 
      //==============================================================
      // Compute X_e from current time till today by solving 
      // the Peebles equation
      //==============================================================

      // The Peebles ODE equation
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };

      // Setting up x-array for the remaining part  
      int npts_peebles = npts_rec_arrays - i + 1;
      Vector x_array_peebles(npts_peebles); 

      for(int j=0; j < npts_peebles; j++){
        x_array_peebles[j] = x_array[i-1+j];
      }

      // Use final Saha value as IC for peebles
      Vector Xe_ic{Xe_arr[i-1]};

      // Solve ODE
      ODESolver peebles_Xe_ode;
      peebles_Xe_ode.solve(dXedx, x_array_peebles, Xe_ic);

      // Fill remainder of arrays 
      auto Xe_arr_peebles = peebles_Xe_ode.get_data_by_component(0);
      for(int j=0; j < npts_peebles; j++){
        Xe_arr[i-1+j] = Xe_arr_peebles[j];
        ne_arr[i-1+j] = Xe_arr_peebles[j] * nb_of_x(x_array_peebles[j]); 
      }

      // Computation finished after completing Peebles 
      break;

    }
  }


  //=============================================================================
  // Splining the result.  
  //=============================================================================
  
  Vector log_Xe_arr = log(Xe_arr);
  Vector log_ne_arr = log(ne_arr);

  log_Xe_of_x_spline.create(x_array, log_Xe_arr, "log_Xe");
  log_ne_of_x_spline.create(x_array, log_ne_arr, "log_ne");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  // Physical constants
  double k_b         = Constants.k_b;
  double G           = Constants.G;
  double m_e         = Constants.m_e;
  double hbar        = Constants.hbar;
  double m_H         = Constants.m_H;
  double epsilon_0   = Constants.epsilon_0;

  // Fetch cosmological parameters
  double Tb          = cosmo->get_TCMB(x);
  double k_b_Tb      = k_b * Tb;
  double nb          = nb_of_x(x); 


  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================

  double pre_factor      = m_e * k_b_Tb / (2.0 * M_PI * hbar * hbar);
  double y               = 1.0/nb * pow(pre_factor, 3.0/2.0) * exp(-epsilon_0 / k_b_Tb);

  // Set Xe=1 at hugehugehuge
  if(y >= 1e9){
    Xe = 1.0;
  } else {
    Xe =  0.5*y * (sqrt(1.0 + 4.0 / y) - 1.0);  
  };

  ne = Xe * nb; 

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of X_e
  double X_e         = Xe[0];
  double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  // const double OmegaB       = cosmo->get_OmegaB(0);
  double Tb           = cosmo->get_TCMB(x);
  double H_of_x       = cosmo->Hp_of_x(x) * exp(-x);
  double nb           = nb_of_x(x); // Assume all baryons are protons 

  // Relevant physical quantities for the Peebles equation
  double hbar_squared     = hbar * hbar;
  double c_squared        = c * c;
  double k_B_Tb           = k_b * Tb;
  double eps0_over_kT     = epsilon_0 / k_B_Tb;
  double e_over_hc_cubed  = epsilon_0*epsilon_0*epsilon_0 / (hbar*hbar*hbar *c*c*c);
  double alpha_squared    = 3.0 * sigma_T * m_e*m_e * c_squared 
                                  / (8.0*M_PI * hbar_squared); 

  //=============================================================================
  // Compute parameters of the Peebles equation 
  //=============================================================================
  double phi_2_of_Tb      = 0.448 * log(eps0_over_kT);
  double alpha_2_of_Tb    = 8.0 / sqrt(3.0 * M_PI) * c * sigma_T
                            * sqrt(eps0_over_kT) * phi_2_of_Tb;
  double beta_of_Tb       = alpha_2_of_Tb 
                            * pow(m_e * k_B_Tb / (2.0*M_PI * hbar_squared), 3.0/2.0)
                            * exp(-eps0_over_kT);
  double nH               = (1.0 - Yp) * nb;
  double n1s              = (1.0 - X_e) * nH;
  double Lambda_alpha     = H_of_x * 27.0 
                            * e_over_hc_cubed
                            / (64.0 * M_PI*M_PI * n1s);

  // Set beta2 to zero to avoid overflow
  // double beta2_of_Tb = 0.0;
  // std::cout << eps0_over_kT << ", ";
  // if(eps0_over_kT < 1.){
    // std::cout << "sug meg: " << x << std::endl;
    // double beta2_of_Tb = beta_of_Tb * exp(3.0 * eps0_over_kT / 4.0);
  // }
  double beta2_of_Tb = exp( - 1.0 * eps0_over_kT / 4.0);

  // std::cout << beta2_of_Tb << ", ";
  
  double C_r_of_Tb        = (lambda_2s1s + Lambda_alpha) 
                            / (lambda_2s1s + Lambda_alpha + beta2_of_Tb); 

  double RHS = C_r_of_Tb/H_of_x * (beta_of_Tb*(1.0 - X_e) - nH*alpha_2_of_Tb * X_e*X_e);

  dXedx[0] = RHS;


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

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  

  // Solve from arbitrary IC
  // Impose tau=0 today manually 
  // Vector tau_ic{1.0};  
  // ODESolver tau_ode;
  // tau_ode.solve(dtaudx, x_array, tau_ic);
  // auto tau_array_forwards = tau_ode.get_data_by_component(0);

  // double tau_today = tau_array_forwards[npts_tau-1]; //  INCORRECT
  // for(int i=0; i<npts_tau; i++){
  //   tau_array_forwards[i] -= tau_today;
  // }

  // Solve backwards, starting from x=0
  Vector tau_ic_back{0.0};   
  ODESolver tau_ode_backwards;
  // tau_ode_backwards.set_accuracy(1e-4, 1e-10, 1e-10);
  tau_ode_backwards.solve(dtaudx, x_backwards_array, tau_ic_back);
  
  // Spline result 
  auto tau_array_reversed = tau_ode_backwards.get_data_by_component(0);
  Vector tau_array(npts_tau);
  for(int i=0; i<npts_tau; i++){
    tau_array[i] = tau_array_reversed[npts_tau-i-1];
  }

  tau_of_x_spline.create(x_array, tau_array, "tau_of_x");
  // tau_of_x_spline.create(x_array, tau_array_forwards, "tau_of_x");



  // Spline derivatives etc. 
  Vector dtau_dx_array(npts_tau);
  Vector g_tilde_array(npts_tau);

  for(int i=0; i<npts_tau; i++){
    dtau_dx_array[i] = - ne_of_x(x_array[i]) * Constants.c * Constants.sigma_T 
                        *exp(x_array[i]) / cosmo->Hp_of_x(x_array[i]);
    g_tilde_array[i] = - dtau_dx_array[i] * exp(-tau_of_x(x_array[i]));

  }

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================


  dtau_dx_spline.create(x_array, dtau_dx_array, "dtau_dx");
  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g_tilde_of_x");

  
  Vector ddtau_ddx_array(npts_tau);
  for(int i=0; i<npts_tau; i++){
    ddtau_ddx_array[i] = dtau_dx_spline.deriv_x(x_array[i]);
  }

  ddtau_ddx_spline.create(x_array, ddtau_ddx_array, "ddtau_ddx");

  //=====================================
  // dg_dx spline may be obsolete. 
  //=====================================
  Vector dg_dx_array(npts_tau);
  
  for(int i=0; i<npts_tau; i++){
    double x_ = x_array[i];
    dg_dx_array[i] = (-ddtauddx_of_x(x_) + dtaudx_of_x(x_)*dtaudx_of_x(x_))
                     * exp(-tau_of_x(x_));
  }

  dg_dx_spline.create(x_array, dg_dx_array, "dg_dx");

  Vector ddg_ddx_array(npts_tau);

  for(int i=0; i<npts_tau; i++){
    // ddtau_ddx_array[i] = dtau_dx_spline.deriv_x(x_array[i]);
    ddg_ddx_array[i] = dg_dx_spline.deriv_x(x_array[i]);
  }

  ddg_ddx_spline.create(x_array, ddg_ddx_array, "ddg_ddx");



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
  // return dg_dx_spline.deriv_x(x);
  return ddg_ddx_spline(x);
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
  return exp(log_Xe_of_x_spline(x));
  // return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
  // return nb_of_x(x) * Xe_of_x(x);
  // return ne_of_x_spline(x);
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
  std::ofstream fp(filename.c_str());

  std::cout << "Writing result to " << filename << std::endl; 

  const int npts       = 9000;
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

