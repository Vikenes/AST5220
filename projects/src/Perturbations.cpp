#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{
  c_  = Constants.c;
  H0_ = cosmo->get_H0();
  H0_squared_ = H0_*H0_;
  OmegaR0_ = cosmo->get_OmegaR();
  OmegaB0_ = cosmo->get_OmegaB();
  OmegaCDM0_ = cosmo->get_OmegaCDM();



}

//====================================================
// Do all the solving
//====================================================
void Perturbations::solve(bool source){

  // Integrate all the perturbation equation and spline the result
  Utils::StartTiming("integrateperturbation");
  integrate_perturbations();
  Utils::EndTiming("integrateperturbation");


  // Compute source functions and spline the result
  if(source){
    Utils::StartTiming("source");
    compute_source_functions();
    Utils::EndTiming("source");
  }

}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================
void Perturbations::integrate_perturbations(){

  //=============================================================
  // Quantities to store 
  //=============================================================
  
  const int nx_nk = n_x * n_k;
  Vector delta_cdm_array_flat(nx_nk);
  Vector delta_b_array_flat(nx_nk);
  Vector v_cdm_array_flat(nx_nk);
  Vector v_b_array_flat(nx_nk);
  Vector Phi_array_flat(nx_nk);
  Vector Psi_array_flat(nx_nk);

  const int N_thetas_store = Constants.n_ell_theta_tc + 1;
  Vector2D Theta_arrays = Vector2D(N_thetas_store, Vector(nx_nk));
  Theta_splines = std::vector<Spline2D>(N_thetas_store);
  const int theta0_idx = 0;//Constants.ind_start_theta_tc;
  const int theta1_idx = 1;//Constants.ind_start_theta_tc + 1;
  const int theta2_idx = 2;//Constants.ind_start_theta_tc + 2;

  //===================================================================
  // Set up the k-array with logarithmic spacing
  //===================================================================
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array(n_k);
  for(int i=0; i<n_k; i++){
    k_array[i] = exp(log_k_array[i]);
  }


  const double hstart = 1e-6;
  const double abserr = 1e-8;
  const double relerr = 1e-8;


  // Loop over all wavenumbers
  #pragma omp parallel for schedule(dynamic, 1)
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    // if( (10*ik) / n_k != (10*ik+10) / n_k ) {
    //   std::cout << (100*ik+100)/n_k << "% " << std::flush;
    //   if(ik == n_k-1) std::cout << std::endl;
    // }

    // Current value of k
    double k = k_array[ik];


    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };
  
    // Find value to integrate to
    auto end_tight  = get_tight_coupling_time(k);
    double x_end_tc = end_tight.first;
    int idx_end_tc  = end_tight.second;

    Vector x_array_tc = Utils::linspace(x_start, x_end_tc, idx_end_tc+1); // CHECK +1 

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // Integrate from x_start -> x_end_tight
    ODESolver ode_tc;
    ode_tc.set_accuracy(hstart, abserr, relerr);
    ode_tc.solve(dydx_tight_coupling, x_array_tc, y_tight_coupling_ini);

    Vector2D y_tc_sol     = ode_tc.get_data();
    Vector2D dydx_tc_sol  = ode_tc.get_derivative_data();
    Vector y_tc_end       = ode_tc.get_final_data();


    //===================================================================
    // Full equation integration
    //===================================================================
   
    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tc_end, x_end_tc, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    int n_x_after_tc = n_x - idx_end_tc;
    Vector x_array_after_tc = Utils::linspace(x_end_tc, x_end, n_x_after_tc); 

    ODESolver ode_after_tc;
    ode_after_tc.set_accuracy(hstart, abserr, relerr);
    ode_after_tc.solve(dydx_full, x_array_after_tc, y_full_ini);

    Vector2D y_after_tc_sol = ode_after_tc.get_data();
    Vector2D dydx_after_tc_sol = ode_after_tc.get_derivative_data();


    for (int ix=0; ix<=idx_end_tc; ix++){
      // Store solution from TC regime 
      double x            = x_array_tc[ix];
      int flat_idx        = ix + n_x * ik;
      Vector y_tc_x       = y_tc_sol[ix];
      Vector dydx_tc_x    = dydx_tc_sol[ix];

      delta_cdm_array_flat[flat_idx]  = y_tc_x[Constants.ind_deltacdm_tc];
      delta_b_array_flat[flat_idx]    = y_tc_x[Constants.ind_deltab_tc];
      v_cdm_array_flat[flat_idx]      = y_tc_x[Constants.ind_vcdm_tc];
      v_b_array_flat[flat_idx]        = y_tc_x[Constants.ind_vb_tc];
      Phi_array_flat[flat_idx]        = y_tc_x[Constants.ind_Phi_tc];

      Theta_arrays[theta0_idx][flat_idx] = y_tc_x[Constants.ind_start_theta_tc];
      Theta_arrays[theta1_idx][flat_idx] = y_tc_x[Constants.ind_start_theta_tc + 1];
      const double Theta2_val               = compute_Theta2_tc(x, k, y_tc_x[Constants.ind_start_theta_tc + 1]);
      Theta_arrays[theta2_idx][flat_idx] = Theta2_val;
      Psi_array_flat[flat_idx]              = compute_Psi(x,k,Theta2_val,y_tc_x[Constants.ind_Phi_tc]); 

    }

    for (int ix=idx_end_tc+1; ix<n_x; ix++){
      // Store solution when TC has ended
      double x                = x_array_full[ix];
      int flat_idx            = ix + n_x * ik;
      int sol_idx             = ix - idx_end_tc;
      Vector y_after_tc_x     = y_after_tc_sol[sol_idx];
      Vector dydx_after_tc_x  = dydx_after_tc_sol[sol_idx];

      delta_cdm_array_flat[flat_idx]  = y_after_tc_x[Constants.ind_deltacdm];
      delta_b_array_flat[flat_idx]    = y_after_tc_x[Constants.ind_deltab];
      v_cdm_array_flat[flat_idx]      = y_after_tc_x[Constants.ind_vcdm];
      v_b_array_flat[flat_idx]        = y_after_tc_x[Constants.ind_vb];
      Phi_array_flat[flat_idx]        = y_after_tc_x[Constants.ind_Phi];

      Theta_arrays[theta0_idx][flat_idx] = y_after_tc_x[Constants.ind_start_theta_tc]; 
      Theta_arrays[theta1_idx][flat_idx] = y_after_tc_x[Constants.ind_start_theta_tc + 1]; 
      Theta_arrays[theta2_idx][flat_idx] = y_after_tc_x[Constants.ind_start_theta_tc + 2]; 

      Psi_array_flat[flat_idx]        = compute_Psi(x, k, y_after_tc_x[Constants.ind_start_theta_tc + 2], Phi_array_flat[flat_idx]);
    }


  }

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  delta_cdm_spline.create(x_array_full, k_array, delta_cdm_array_flat);
  delta_b_spline.create(x_array_full, k_array, delta_b_array_flat);
  v_cdm_spline.create(x_array_full, k_array, v_cdm_array_flat);
  v_b_spline.create(x_array_full, k_array, v_b_array_flat);
  Phi_spline.create(x_array_full, k_array, Phi_array_flat);
  Psi_spline.create(x_array_full, k_array, Psi_array_flat);

  for(int ell=0; ell<Theta_splines.size(); ell++){
    std::string name = "theta" + std::to_string(ell) + "spline"; 
    Theta_splines[ell].create(x_array_full, k_array, Theta_arrays[ell], name);
  }


}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  // I THINK THESE ARE OBSOLETE  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];

  // Scalar quantities (Gravitational potential, baryons and CDM)
  double Psi  = -2.0 / 3.0;
  Phi         = - Psi;
  delta_cdm   = -1.5 * Psi;
  delta_b     = -1.5 * Psi; 
  v_cdm       = - c_ * k / (2.0 * cosmo->Hp_of_x(x)) * Psi;
  v_b         = v_cdm;    

  // Photon temperature perturbations (Theta_ell)
  Theta[0]    = -0.5 * Psi;
  Theta[1]    = - v_cdm / 3.0;


  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================
Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double &Phi             =  y[Constants.ind_Phi];
  double *Theta           = &y[Constants.ind_start_theta];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  const double ck_over_Hp_dtaudx = c_ * k / (cosmo->Hp_of_x(x) * rec->dtaudx_of_x(x));

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  Phi       = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b   = delta_b_tc;
  v_cdm     = v_cdm_tc;
  v_b       = v_b_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  Theta[2] = compute_Theta2_tc(x, k, Theta[1]);
  for (int l=3; l<n_ell_theta; l++){
    Theta[l] = -(double)l / (double)(2.*l+1.) * ck_over_Hp_dtaudx * Theta[l-1];
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================
std::pair<double,int> Perturbations::get_tight_coupling_time(const double k) const{
  // Ensure fuckup if end time is not found 
  double x_tc_end = 1000.0;
  int  idx_tc_end = n_x*10;

  //=============================================================================
  // Remember all the three conditions in Callin
  //=============================================================================
  for (int i=0; i<n_x; i++){
    double x          = x_array_full[i];
    double ck         = c_ * k;
    double tau_prime  = rec->dtaudx_of_x(x);
    tau_prime         *= -1; // Make dtau/dx a positive quantity 

    bool cond1 = tau_prime < 10.0;
    bool cond2 = tau_prime < 10.0 * ck / cosmo->Hp_of_x(x);
    bool cond3 = x > -8.3;  
    // bool cond3 = rec->Xe_of_x(x) < 0.99;  

    if (cond1 || cond2 || cond3) {
      // Tight-coupling regime not valid 
      idx_tc_end = i - 1;
      x_tc_end = x_array_full[idx_tc_end];
      break; 
    }
  }

  return std::pair<double,int>(x_tc_end, idx_tc_end);
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================

  // Use quadratically distributed k-values
  const double kmin = Constants.k_min;
  const double kmax = Constants.k_max;

  Vector log_k_array = Utils::linspace(log(kmin), log(kmax), n_k);
  Vector k_array = exp(log_k_array); 

  // for(int i=0; i<n_k; i++){
    // double i_term = (double)i / (double)(n_k - 1.0);
    // k_array[i] = kmin + Delta_k * i_term*i_term;
  // }  
  Vector x_array = x_array_full;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){

    const double x = x_array[ix];

    // Fetch cosmo quantities  
    const double Hp           = cosmo->Hp_of_x(x);
    const double dHp_dx       = cosmo->dHpdx_of_x(x);
    const double ddHp_ddx     = cosmo->ddHpddx_of_x(x);


    // Fetch rec quantities  
    const double tau          = rec->tau_of_x(x);
    const double g_tilde      = rec->g_tilde_of_x(x);
    const double dgdx_tilde   = rec->dgdx_tilde_of_x(x);
    const double ddgddx_tilde = rec->ddgddx_tilde_of_x(x);

    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;


      // Constants
      const double ck           = c_ * k; 

      // Fetch pert quantities 
      const double Theta0       = get_Theta(x, k, 0);
      const double Pi           = get_Theta(x, k, 2);
      const double Psi          = get_Psi(x,k);
      const double v_b          = get_v_b(x,k);

      const double dPsi_dx      = Psi_spline.deriv_x(x,k);
      const double dPhi_dx      = Phi_spline.deriv_x(x,k);
      const double dvb_dx       = v_b_spline.deriv_x(x,k);
      const double dPi_dx       = Theta_splines[2].deriv_x(x,k);
      const double ddPi_ddx     = Theta_splines[2].deriv_xx(x,k);

      // Compute terms in Source function. 
      const double first_term         = g_tilde * (Theta0 + Psi + Pi/4.0);
      const double second_term        = exp(-tau) * (dPsi_dx - dPhi_dx);
      const double d_Hp_gtilde_vb_dx  = dHp_dx * g_tilde * v_b 
                                        + Hp * dgdx_tilde * v_b 
                                        + Hp * g_tilde * dvb_dx;
      const double third_term         = - d_Hp_gtilde_vb_dx / ck;

      const double fourth_term_1      = g_tilde * Pi * (dHp_dx*dHp_dx + Hp * ddHp_ddx);
      const double fourth_term_2      = 3.0 * Hp * dHp_dx * (dgdx_tilde*Pi + g_tilde * dPi_dx);
      const double fourth_term_3      = Hp*Hp * (ddgddx_tilde*Pi + 2.0 * dgdx_tilde * dPi_dx + g_tilde*ddPi_ddx);
      const double fourth_term        = 3.0 / (4.0 * ck*ck) * (fourth_term_1+fourth_term_2+fourth_term_3);  


      // Temperatur source
      ST_array[index] = first_term + second_term + third_term + fourth_term;
    }
  }

  // Spline the source functions
  ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");
  

}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];


  const double a          = exp(x);
  const double ck         = c_ * k;
  // const double H0         = cosmo->get_H0();
  const double Hp         = cosmo->Hp_of_x(x);

  const double a_squared  = a*a;
  const double ck_squared = ck * ck; 
  // const double H0_squared = H0 * H0;
  const double Hp_squared = Hp * Hp;

  const double dHp_dx     = cosmo->dHpdx_of_x(x);
  const double ck_over_Hp = ck / Hp;
  const double dtau_dx    = rec->dtaudx_of_x(x);
  const double ddtau_ddx  = rec->ddtauddx_of_x(x);

  // const double OmegaCDM0  = cosmo->get_OmegaCDM();
  // const double OmegaB0    = cosmo->get_OmegaB();
  // const double OmegaR0    = cosmo->get_OmegaR();
  const double R          = 4.0 * OmegaR0_ / (3.0 * OmegaB0_ * a);



  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================
  const double Theta2 = compute_Theta2_tc(x, k, Theta[1]);
  const double Psi    = compute_Psi(x,k,Theta2,Phi);

  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx = Psi - ck_over_Hp*ck_over_Hp * Phi / 3.0 
              + H0_squared_ / (2.0 * Hp_squared)
              * ((OmegaCDM0_ * delta_cdm + OmegaB0_ * delta_b) / a + 4.0 * OmegaR0_ * Theta[0] / a_squared);

  ddelta_cdmdx = ck_over_Hp * v_cdm - 3.0 * dPhidx;
  ddelta_bdx   = ck_over_Hp * v_b   - 3.0 * dPhidx;
  dv_cdmdx     = -v_cdm - ck_over_Hp * Psi;

  // SET: Photon multipoles (Theta_ell)

  dThetadx[0] = -ck_over_Hp * Theta[1] - dPhidx;
  const double q_nominator = - ((1.0 - R)*dtau_dx + (1.0 + R)*ddtau_ddx)*(3.0*Theta[1] + v_b)
                             - ck_over_Hp*Psi 
                             + (1.0 - dHp_dx/Hp) * ck_over_Hp * (-Theta[0] + 2.0*Theta2) 
                             - ck_over_Hp * dThetadx[0];
  const double q_denominator = (1.0 + R)*dtau_dx + dHp_dx/Hp - 1.0;
  

  const double q  = q_nominator / q_denominator;
  dv_bdx          = (-v_b-ck_over_Hp*Psi + R*(q + ck_over_Hp*(-Theta[0] + 2.0*Theta2) - ck_over_Hp*Psi )) / (1.0 + R);
  dThetadx[1]     = (q - dv_bdx) / 3.0;


  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================
int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];

  // Cosmological parameters and variables
  // const double H0             = cosmo->get_H0();
  // const double H0_squared     = H0 * H0;
  const double Hp             = cosmo->Hp_of_x(x);
  const double a              = exp(x);
  const double a_inv          = 1.0 / a; 
  const double a_inv_squared  = a_inv * a_inv;
  // const double OmegaR0        = cosmo->get_OmegaR();
  // const double OmegaCDM0      = cosmo->get_OmegaCDM();
  // const double OmegaB0        = cosmo->get_OmegaB();
  const double eta            = cosmo->eta_of_x(x);

  // Recombination variables
  const double dtau_dx = rec->dtaudx_of_x(x);

  // Parameters 
  // const double c            = Constants.c;
  const double ck           = c_ * k;
  const double ck_over_Hp   = ck / Hp;
  const double ck_over_3Hp  = ck_over_Hp / 3.0; 
  const double R            = 4.0 * OmegaR0_ * a_inv / (3.0 * OmegaB0_);

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  const double Psi = compute_Psi(x,k,Theta[2],Phi);


  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx = Psi - ck_over_Hp*ck_over_Hp * Phi / 3.0 
                + H0_squared_ / (2.0 * Hp*Hp) *
                (OmegaCDM0_ * delta_cdm * a_inv 
                 + OmegaB0_ * delta_b * a_inv
                 + 4.0 * OmegaR0_ * Theta[0] * a_inv_squared);

  ddelta_cdmdx  = ck_over_Hp * v_cdm - 3.0 * dPhidx;
  ddelta_bdx    = ck_over_Hp * v_b - 3.0 * dPhidx;
  dv_cdmdx      = -v_cdm - ck_over_Hp * Psi;
  dv_bdx        = -v_b - ck_over_Hp * Psi + dtau_dx * R * (3.0 * Theta[1] + v_b);

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -ck_over_Hp * Theta[1] - dPhidx;
  dThetadx[1] = ck_over_3Hp * Theta[0] - 2.0 * ck_over_3Hp * Theta[2] + ck_over_3Hp * Psi 
                + dtau_dx * (Theta[1] + v_b/3.0);
                
  const int l_max = n_ell_theta - 1;

  for (int l=2; l<l_max; l++){
    double ell_denom = (2.*l + 1.);
    // double ell_minus_term = l * ck_over_Hp * Theta[l-1] / ell_denom;
    // double ell_plus_term  = (l + 1.) * ck_over_Hp * Theta[l+1] / ell_denom;
    // double ell_term = dtau_dx * Theta[l];
    dThetadx[l] = l * ck_over_Hp * Theta[l-1] / ell_denom 
                  - (l + 1.) * ck_over_Hp * Theta[l+1] / ell_denom 
                  + dtau_dx * Theta[l];
    if(l==2){
      dThetadx[l] -= dtau_dx * Theta[l]/10.0;
    } 

  }
  dThetadx[l_max] = ck_over_Hp * Theta[l_max - 1]
                              - c_ * (l_max + 1.0) / (Hp * eta) * Theta[l_max]
                              + dtau_dx * Theta[l_max];


  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return get_Theta(x,k, 2);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_splines[ell](x,k);
}
double Perturbations::compute_Theta2_tc(const double x, const double k, const double Theta1) const{
  return - 20.0 * c_*k / (45.0 * cosmo->Hp_of_x(x) * rec->dtaudx_of_x(x)) * Theta1;
}
double Perturbations::compute_Psi(const double x, const double k, const double Theta2, const double Phi) const{
  // const double H0 = cosmo->get_H0();
  const double ck = c_ * k;
  return -Phi - 12.0 * H0_squared_ / (exp(2.0*x) * ck*ck) * OmegaR0_ * Theta2;
}

//====================================================
// Print some useful info about the class
//====================================================
void Perturbations::info() const{
  /*
    A hot mess at the moment. Clean up later
  */
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================
void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_v_b(x,k)       << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };

  std::cout << "Saving n=" << npts << " data points to '" << filename << "'" << std::endl;

  double x_mr_eq = cosmo->get_mr_equality(x_array);
  auto end = get_tight_coupling_time(k);
  double xend = end.first; 
  double x_entry = cosmo->get_k_eta_equals_unity(k);
  fp << xend << " " << x_mr_eq << " " << x_entry << "\n";
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

