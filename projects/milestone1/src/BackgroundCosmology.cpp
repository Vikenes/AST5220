#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0, ...
  //=============================================================================
  H0 = Constants.H0_over_h * h; 
  double kT = Constants.k_b * TCMB;

  double hbar = Constants.hbar;
  double c = Constants.c;
  double hbar3_c5 = hbar*hbar*hbar * c*c*c*c*c;

  OmegaR = M_PI*M_PI / 15.0 * kT*kT*kT*kT / hbar3_c5 * 8.0 * M_PI * Constants.G / (3.0 * H0*H0);
  
  OmegaNu = Neff * 7.0 / 8.0 * pow(4.0/11.0, 4.0/3.0) * OmegaR;
  OmegaLambda = 1 - (OmegaB + OmegaCDM + OmegaK + OmegaR + OmegaNu); 

  

}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  // Utils::StartTiming("Eta");
  //=============================================================================
  // Solve the ODE for eta(x) and t(x), and create a Spline. 
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, nx);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    detadx[0] = Constants.c / Hp_of_x(x);
    return GSL_SUCCESS;
  };

  // The initial conditions for deta/dx
  double eta_ini = Constants.c / Hp_of_x(x_start);
  Vector eta_ic{eta_ini};

  // Solve the ODE 
  ODESolver eta_ode;
  eta_ode.solve(detadx, x_array, eta_ic);

  // Spline result 
  auto eta_array = eta_ode.get_data_by_component(0);
  eta_of_x_spline.create(x_array, eta_array, "Eta(x) spline");

  // Utils::EndTiming("Eta");


  // Utils::StartTiming("t");
  //=============================================================================
  // Solve ODE for t(x) and Spline result.
  //=============================================================================

  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
    dtdx[0] = 1.0 / H_of_x(x);
    return GSL_SUCCESS;
  };

  double t_ini = 1.0 / (2.0 * H_of_x(x_start));
  Vector t_ic{t_ini};

  // Solve ODE 
  ODESolver t_ode;
  t_ode.solve(dtdx, x_array, t_ic);

  // Spline 
  auto t_array = t_ode.get_data_by_component(0);
  t_of_x_spline.create(x_array, t_array, "t(x) spline");

  // Utils::EndTiming("t");
}


// Solve the background
void BackgroundCosmology::solve_eta(){
  // Utils::StartTiming("Eta");
  //=============================================================================
  // Solve the ODE for eta(x), and create a Spline. 
  // Use for supernova fitting, to avoid unnecessary computations of t(x).
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, nx);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    detadx[0] = Constants.c / Hp_of_x(x);
    return GSL_SUCCESS;
  };

  // The initial conditions for deta/dx
  double eta_ini = Constants.c / Hp_of_x(x_start);
  Vector eta_ic{eta_ini};

  // Solve the ODE 
  ODESolver eta_ode;
  eta_ode.solve(detadx, x_array, eta_ic);

  // Spline result 
  auto eta_array = eta_ode.get_data_by_component(0);
  eta_of_x_spline.create(x_array, eta_array, "Eta(x) spline");

  // Utils::EndTiming("Eta");


}

//====================================================
// Get methods
//====================================================

/*Computes H(x)*/
double BackgroundCosmology::H_of_x(double x) const{

  double Hx = H0 * sqrt((OmegaB + OmegaCDM)*exp(-3*x) 
                        + (OmegaR + OmegaNu) * exp(-4*x)
                        + OmegaK * exp(-2*x) 
                        + OmegaLambda);
  return Hx;
}

/*Computes Hp(x)*/
double BackgroundCosmology::Hp_of_x(double x) const{
  // Reduce number of exp-calls.
  // Increase efficiency for supernova fitting. 
  double exp_of_minus_x = exp(-x); 
  double exp_of_minus_2x = exp_of_minus_x * exp_of_minus_x;

  double Hp = H0 * sqrt((OmegaB + OmegaCDM) * exp_of_minus_x 
                        + (OmegaR + OmegaNu) * exp_of_minus_2x 
                        + OmegaK 
                        + OmegaLambda / exp_of_minus_2x);
  return Hp;
}

/*
  Compute Hp'(x). Using that Hp'(x)=1/(2Hp)*d/dx(...)
*/
double BackgroundCosmology::dHpdx_of_x(double x) const{
  // Term inside square root. 
  double dv_sqrt_term = -      (OmegaB + OmegaCDM) * exp(-x)
                        - 2 * (OmegaR + OmegaNu) * exp(-2*x) 
                        + 2 * OmegaLambda * exp(2*x);

  double dHp = pow(H0,2) / (2 * Hp_of_x(x)) * dv_sqrt_term; 

  return dHp;
}

/*Compute Hp''(x). Using chain rule on Hp'(x) and simplifying.*/
double BackgroundCosmology::ddHpddx_of_x(double x) const{
  // Derivative of square root term 
  double dv2_sqrt_term = (OmegaB + OmegaCDM) * exp(-x)
                        + 4 * (OmegaR + OmegaNu) * exp(-2*x) 
                        + 4 * OmegaLambda * exp(2*x);

  // Derivative of 1/Hp, multiplied by square root term and simplified.
  double dHp_term = - 2 / pow(H0,2) * pow(dHpdx_of_x(x),2);

  double ddHp = pow(H0,2) / (2 * Hp_of_x(x))
                * (dv2_sqrt_term + dHp_term);                       

  return ddHp;
}

/*OmegaB(x)*/
double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;
  return OmegaB / (exp(x) * pow(Hp_of_x(x) / H0 , 2) );
}

/*OmegaR(x)*/
double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;
  return OmegaR / (exp(2*x) * pow(Hp_of_x(x) / H0 , 2) );
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  return OmegaNu / (exp(2*x) * pow(Hp_of_x(x) / H0 , 2) );
}

/*OmegaCDM(x)*/
double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;
  return OmegaCDM / (exp(x) * pow(Hp_of_x(x) / H0 , 2) );
}

/*OmegaLambda(x)*/
double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;
  return OmegaLambda / pow(H_of_x(x) / H0,2);
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;
  return OmegaK / pow(Hp_of_x(x) / H0, 2);
}


//=============================================================================
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  double dL = get_comoving_distance_of_x(x) * exp(-x);
  return dL;
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  double chi = eta_of_x(0.0) - eta_of_x(x);
  if(OmegaK == 0.0) return chi;

  double OmegaK_factor = sqrt(fabs(OmegaK)) * H0 * chi / Constants.c; 
  if(OmegaK < 0) return chi * sin(OmegaK_factor) / OmegaK_factor;
  if(OmegaK > 0) return chi * std::sinh(OmegaK_factor) / OmegaK_factor;

  return 0.0;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = x_start;
  const double x_max = x_end;
  const int    n_pts = npts;

  std::cout << "Writing to " << filename << ", for " 
  << x_min << " < x < " << x_max << std::endl;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";
    fp << eta_of_x(x)        << " ";
    fp << get_t_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << get_luminosity_distance_of_x(x) << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

