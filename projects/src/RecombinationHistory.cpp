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
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;
      

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;

      // std::cout << "Xe=" << Xe_arr[i] <<", x=" << x_array[i] << std::endl;
      // if(x_array[i]>-4.3) break;

    } else {
      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      std::cout << "Peebles entered" << std::endl;
      // The Peebles ODE equation
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      //==============================================*c===============================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================

      // create x-array for remaining region
      int npts_peebles = npts_rec_arrays - i;
      std::cout << npts_peebles << std::endl;
      Vector x_array_peebles(npts_peebles); 
      for(int j=0; j < npts_peebles; j++){
        x_array_peebles[j] = x_array[i+j];
      }

      // Use final Saha result as IC for peebles
      Vector Xe_ic{Xe_arr[i-1]};

      // Solve ODE
      ODESolver peebles_Xe_ode;
      peebles_Xe_ode.solve(dXedx, x_array_peebles, Xe_ic);

      // Fill remainder of arrays 
      auto Xe_arr_peebles = peebles_Xe_ode.get_data_by_component(0);
      for(int j=0; j < npts_peebles; j++){
        Xe_arr[i+j] = Xe_arr_peebles[j];
        ne_arr[i+j] = Xe_arr_peebles[j] * nb_of_x(x_array_peebles[j]); 
      }

      // Computation finished after completing Peebles 
      break;

    }
  }


  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  
  Vector log_Xe_arr = log(Xe_arr);
  Vector log_ne_arr = log(ne_arr);

  // std::cout << "Length of     x_ = " << x_array.size() << std::endl;
  // std::cout << "Length of log Xe = " << log_Xe_arr.size() << std::endl;
  // std::cout << "Length of     Xe = " << Xe_arr.size() << std::endl;


  log_Xe_of_x_spline.create(x_array, log_Xe_arr, "log_Xe");
  log_ne_of_x_spline.create(x_array, log_ne_arr, "log_ne");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  // double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  // const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  // const double OmegaB      = cosmo->get_OmegaB(0);
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
  double y               = 1.0/nb * pow(pre_factor, 1.5) * exp(-epsilon_0 / k_b_Tb);

  // Set Xe=1 at hugehugehuge
  if(y > 1e7){
    Xe = 1.0;
  } else {
    // Xe =  pow(y*y + 4.0 * y, 0.5)/2.0 ;  
    Xe = (pow(1 + 4.0 / y, 0.5) - 1) * y/2.0;  
  };

  ne = Xe * nb; 

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

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
  const double OmegaB       = cosmo->get_OmegaB(0);
  double Tb                 = cosmo->get_TCMB(x);
  double H_of_x             = cosmo->H_of_x(x);
  const double H0           = cosmo->get_H0();

  // Relevant physical quantities for the Peebles equation
  const double hbar_squared     = hbar * hbar;
  const double c_squared        = c * c;
  double k_B_Tb           = k_b * Tb;
  const double eps0_over_kT     = epsilon_0 / k_B_Tb;
  const double e_over_hc_cubed  = epsilon_0*epsilon_0*epsilon_0 / (hbar*hbar*hbar *c*c*c);
  const double alpha_squared    = 3.0 * sigma_T * m_e*m_e * c_squared 
                                  / (8.0*M_PI * hbar_squared); 

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  double phi_2_of_Tb      = 0.448 * log(eps0_over_kT);
  double alpha_2_of_Tb    = 64.0 * M_PI / (sqrt(27.0 * M_PI))
                            * alpha_squared * hbar_squared / (m_e*m_e*c)
                            * sqrt(eps0_over_kT) * phi_2_of_Tb;
  double beta_of_Tb       = alpha_2_of_Tb 
                            * pow(m_e * k_B_Tb / (2.0*M_PI * hbar_squared), 3.0/2.0)
                            * exp(-eps0_over_kT);
  double nH               = 3.0 * H0*H0 * OmegaB / (8.0*M_PI * G * m_H * a*a*a);
  double n1s              = (1.0 - X_e) * nH;
  double Lambda_alpha     = H_of_x * 27.0 
                            * e_over_hc_cubed
                            / (64.0 * M_PI*M_PI * n1s);

  double beta2_of_Tb = 0.0;
  if(eps0_over_kT < 200){
    double beta2_of_Tb = beta_of_Tb * exp(3.0 * eps0_over_kT / 4.0);
  }
  
  double C_r_of_Tb        = (lambda_2s1s + Lambda_alpha) 
                            / (lambda_2s1s + Lambda_alpha + beta2_of_Tb); 

  double RHS = C_r_of_Tb/H_of_x * (beta_of_Tb*(1.0 - X_e) - nH*alpha_2_of_Tb * X_e*X_e);

  // dXedx[0] = 0.0;
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
  Vector x_array = Utils::linspace(x_start, x_end, npts_tau);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    const double c       = Constants.c; 
    const double sigma_T = Constants.sigma_T;
    dtaudx[0] = - c * sigma_T * ne_of_x(x) / cosmo->H_of_x(x); // Arbitrary IC
    // dtaudx[0] = c * sigma_T * ne_of_x(x) / cosmo->H_of_x(x);   // For x=0 as IC
    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  
  // Attempting x=0 as IC 
  // Vector tau_ic{0.0};

  // Arbitrary IC
  Vector tau_ic{1.0};
  
  ODESolver tau_ode;
  tau_ode.solve(dtaudx, x_array, tau_ic);

  // Spline result 


  // Arbitrary IC
  auto tau_array = tau_ode.get_data_by_component(0);


  double tau_today = tau_array[npts_tau-1];  
  for(int i=0; i<npts_tau; i++){
    tau_array[i] -= tau_today;
  }


  tau_of_x_spline.create(x_array, tau_array, "tau_of_x");


  // Spline derivatives etc. 
  Vector dtau_dx_array(npts_tau);
  Vector g_tilde_array(npts_tau);

  for(int i=0; i<npts_tau; i++){
    dtau_dx_array[i] = - ne_of_x(x_array[i]) * Constants.c * Constants.sigma_T 
                        / cosmo->H_of_x(x_array[i]);
    g_tilde_array[i] = - dtau_dx_array[i] * exp(-tau_of_x(x_array[i]));

  }
  // Vector ddtau_dxx_array(npts); 
  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================

  // for(int i=0; i<npts_tau; i++){
    // g_tilde_array[i] = - dtaudx_of_x(x_array[i]) * exp(-tau_of_x(x_array[i]));
  // } 

  dtau_dx_spline.create(x_array, dtau_dx_array, "dtau_dx");
  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g_tilde_of_x");


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

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  // 
  // Analytical value 
  // double dtau = - ne_of_x(x) * Constants.c * Constants.sigma_T / cosmo->H_of_x(x);
  //
  // return tau_of_x_spline.deriv_x(x);
  return dtau_dx_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return dtau_dx_spline.deriv_x(x);;
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


double RecombinationHistory::nb_of_x(double x) const{
  //=============================================================================
  // Compute baryon number density at a given x
  //=============================================================================
  double a_cubed_inv = exp(-3.0*x);
  double OmegaB      = cosmo->get_OmegaB(0);  
  // double nb          = 

  return OmegaB * cosmo->get_rho_c0() * a_cubed_inv / Constants.m_H;;
}


double RecombinationHistory::Xe_of_x(double x) const{
  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
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

  const int npts       = nx;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  fp << "nx: " << npts << ", npts_rec: " << npts_rec_arrays << ", npts_tau_g: " << npts_tau << "\n";
  fp << "x Xe ne tau dtaudx ddtauddx g dgdx ddgddx \n";
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
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

