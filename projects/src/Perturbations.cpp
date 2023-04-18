#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================
void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================
void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //=============================================================
  // Quantities to store 
  //=============================================================
  Vector nk_vector(n_k);
  Vector2D delta_cdm_array(n_x, nk_vector);
  Vector2D delta_b_array(n_x, nk_vector);
  Vector2D v_cdm_array(n_x, nk_vector);
  Vector2D v_b_array(n_x, nk_vector);
  Vector2D Phi_array(n_x, nk_vector);
  Vector2D Psi_array(n_x, nk_vector);

  Vector2D Theta0_array(n_x, nk_vector);
  Vector2D Theta1_array(n_x, nk_vector);
  Vector2D Theta2_array(n_x, nk_vector);


  Vector delta_cdm_array_flat(n_x * n_k);
  Vector delta_b_array_flat(n_x * n_k);
  Vector v_cdm_array_flat(n_x * n_k);
  Vector v_b_array_flat(n_x * n_k);
  Vector Phi_array_flat(n_x * n_k);
  Vector Psi_array_flat(n_x * n_k);

  Vector Theta0_array_flat(n_x * n_k);
  Vector Theta1_array_flat(n_x * n_k);
  Vector Theta2_array_flat(n_x * n_k);



  //===================================================================
  // Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with a
  // logarithmic spacing
  //===================================================================
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array(n_k);
  for(int i=0; i<n_k; i++){
    k_array[i] = exp(log_k_array[i]);
  }


  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    auto end_tight = get_tight_coupling_time(k);
    double x_end_tc = end_tight.first;
    int idx_end_tc  = end_tight.second;

    Vector x_array_tc = Utils::linspace(x_start, x_end_tc, idx_end_tc + 1); // CHECK +1 

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    ODESolver ode_tc;
    ode_tc.solve(dydx_tight_coupling, x_array_tc, y_tight_coupling_ini);

    Vector2D y_tc_sol     = ode_tc.get_data();
    Vector2D dydx_tc_sol  = ode_tc.get_derivative_data();
    Vector y_tc_end       = ode_tc.get_final_data();

    double ck = Constants.c * k;
    double ck_squared = ck * ck;
    const double H0 = cosmo->get_H0();
    const double H0_squared = H0 * H0;
    const double OmegaR0 = cosmo->get_OmegaR();

    for (int ix=0; ix<idx_end_tc; ix++){
      double x = x_array_tc[ix];
      double a_squared = exp(2*x);
      double Hp       = cosmo->Hp_of_x(x);
      double dtau_dx  = rec->dtaudx_of_x(x);

      delta_cdm_array[ix][ik] = y_tc_sol[ix][Constants.ind_deltacdm_tc];
      delta_b_array[ix][ik]   = y_tc_sol[ix][Constants.ind_deltab_tc];
      v_cdm_array[ix][ik]     = y_tc_sol[ix][Constants.ind_vcdm_tc];
      v_b_array[ix][ik]       = y_tc_sol[ix][Constants.ind_vb_tc];
      Theta0_array[ix][ik]    = y_tc_sol[ix][Constants.ind_start_theta_tc];

      double Phi_tc           = y_tc_sol[ix][Constants.ind_Phi_tc];
      double Theta1_tc        = y_tc_sol[ix][Constants.ind_start_theta_tc + 1];
      double Theta2_tc        = - 20.0 * ck * Theta1_tc / (45.0 * Hp * dtau_dx);
      Phi_array[ix][ik]       = Phi_tc;
      Theta1_array[ix][ik]    = Theta1_tc;
      Theta2_array[ix][ik]    = Theta2_tc;
      Psi_array[ix][ik]       = - Phi_tc - 12.0*H0_squared/(ck_squared*a_squared)*OmegaR0 * Theta2_tc;
    }



    //===================================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
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
    ode_after_tc.solve(dydx_full, x_array_after_tc, y_full_ini);

    Vector2D y_after_tc_sol = ode_after_tc.get_data();
    Vector2D dydx_after_tc_sol = ode_after_tc.get_derivative_data();

    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //
    //===================================================================
    //...
    //...

  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...
  // ...
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
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
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================
  // ...
  // ...

  
  // Scalar quantities (Gravitational potential, baryons and CDM)
  double Psi  = -2.0 / 3.0;
  Phi         = - Psi;
  delta_cdm   = -1.5 * Psi;
  delta_b     = -1.5 * Psi; 
  v_cdm       = - Constants.c * k / (2.0 * cosmo->Hp_of_x(x)) * Psi;
  v_b         = v_cdm;    

  // Photon temperature perturbations (Theta_ell)
  Theta[0]    = -0.5 * Psi;
  Theta[1]    = v_cdm / 3.0;


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
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  double ck_over_Hp_dtaudx = Constants.c * k / (cosmo->Hp_of_x(x) * rec->dtaudx_of_x(x));

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  Phi       = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b   = delta_b_tc;
  v_cdm     = v_cdm_tc;
  v_b       = v_b_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  Theta[2] = -4.0 / 9.0 * ck_over_Hp_dtaudx * Theta[1];
  for (int l=3; l<n_ell_theta; l++){
    Theta[l] = -(double)l / (2.*l+1.) * ck_over_Hp_dtaudx * Theta[l-1];
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================
std::pair<double,int> Perturbations::get_tight_coupling_time(const double k) const{
  double x_tight_coupling_end = 0.0;
  int  idx_tight_coupling_end = 0;

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  for (int i=0; i<n_x; i++){
    double x = x_array_full[i];
    double tau_prime = rec->dtaudx_of_x(x);
    tau_prime *= -1; // Make dtau/dx<0 a positive quantity 
    if (tau_prime < 10 || tau_prime < 10 * Constants.c * k / cosmo->Hp_of_x(x) || x > -8.3) {
      // Tight-coupling regime not valid 

      // DOUBLE CHECK THAT IT'S NOT X-DX AND I-1
      x_tight_coupling_end = x; 
      idx_tight_coupling_end = i;
      break; 
    }
  }

  return std::pair<double,int>(x_tight_coupling_end, idx_tight_coupling_end);
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  Vector k_array;
  Vector x_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);
      // ...
      // ...

      // Temperatur source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
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
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];


  const double a          = exp(x);
  const double ck         = Constants.c * k;
  const double H0         = cosmo->get_H0();
  const double Hp         = cosmo->Hp_of_x(x);

  const double a_squared  = a*a;
  const double ck_squared = ck * ck; 
  const double H0_squared = H0 * H0;
  const double Hp_squared = Hp * Hp;

  const double dHp_dx     = cosmo->dHpdx_of_x(x);
  const double ck_over_Hp = ck / Hp;
  const double dtau_dx    = rec->dtaudx_of_x(x);
  const double ddtau_ddx  = rec->ddtauddx_of_x(x);

  const double OmegaCDM0  = cosmo->get_OmegaCDM();
  const double OmegaB0    = cosmo->get_OmegaB();
  const double OmegaB0    = cosmo->get_OmegaB();
  const double OmegaR0    = cosmo->get_OmegaR();
  const double R          = 4.0 * OmegaR0 / (3.0 * OmegaB0 * a);




  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================
  const double Theta2 = - 20.0 * ck_over_Hp / (45.0 * dtau_dx) * Theta[1];
  const double Psi    = - Phi - 12.0 * H0_squared * OmegaR0 * Theta2 / (ck_squared * a_squared);
  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx = Psi - ck_over_Hp*ck_over_Hp * Phi / 3.0 
              + H0_squared / (2.0 * Hp_squared)
              * ((OmegaCDM0 * delta_cdm + OmegaB0 * delta_b) / a + 4.0 * OmegaR0 * Theta[0] / a_squared);

  ddelta_cdmdx = ck_over_Hp * v_cdm - 3.0 * dPhidx;
  ddelta_bdx   = ck_over_Hp * v_b   - 3.0 * dPhidx;
  dv_cdmdx     = -v_cdm - ck_over_Hp * Psi;

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -ck_over_Hp * Theta[1] - dPhidx;

  const double q_nominator = - ((1.0 - R)*dtau_dx + (1.0 + R)*ddtau_ddx)*(3.0*Theta[1] + v_b)
                             - ck_over_Hp*Psi 
                             + (1 - dHp_dx/Hp) * ck_over_Hp * (-Theta[0] + 2.0*Theta2) 
                             - ck_over_Hp * dThetadx[0];
  const double q_denominator = (1.0 + R)*dtau_dx + dHp_dx/Hp - 1;
  

  const double q  = q_nominator / q_denominator;
  dv_bdx          = (-v_b - ck_over_Hp * Psi + R*(q + ck_over_Hp*(-Theta[0] + 2.0*Theta2) - ck_over_Hp*Psi )) / (1.0 + R);
  dThetadx[1]     = (q - dv_bdx) / 3.0;


  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================
int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables
  const double H0 = cosmo->get_H0();
  const double H0_squared = H0 * H0;
  const double Hp = cosmo->Hp_of_x(x);
  const double a_inv = exp(-x);
  const double a_inv_squared = a_inv * a_inv;
  const double OmegaR0 = cosmo->get_OmegaR();
  const double OmegaCDM0 = cosmo->get_OmegaCDM();
  const double OmegaB0   = cosmo->get_OmegaB();
  const double eta = cosmo->eta_of_x(x);
  

  // Recombination variables
  const double dtau_dx = rec->dtaudx_of_x(x);

  // Parameters 
  const double c = Constants.c;
  const double ck = c * k;
  const double ck_over_Hp = ck / Hp;
  const double ck_over_3Hp = ck_over_Hp / 3.0; 
  const double R = 4.0 * OmegaR0 * a_inv / (3.0 * OmegaB0);

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  const double Psi = -Phi - 12.0*H0_squared*OmegaR0*Theta[2]*a_inv_squared / (ck*ck);

  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx = Psi - ck_over_Hp*ck_over_Hp * Phi / 3.0 
                + H0_squared / (2.0 * Hp*Hp) *
                (OmegaCDM0 * delta_cdm * a_inv 
                 + OmegaB0 * delta_b * a_inv
                 + 4.0 * OmegaR0 * Theta[0] * a_inv_squared);

  ddelta_cdmdx = ck_over_Hp * v_cdm - 3.0 * dPhidx;
  ddelta_bdx = ck_over_Hp * v_b - 3.0 * dPhidx;
  dv_cdmdx = -v_cdm - ck_over_Hp * Psi;
  dv_bdx = -v_b - ck_over_Hp * Psi + dtau_dx * R * (3.0 * Theta[1] + v_b);

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
                              - c * (l_max + 1.0) / (Hp * eta) * Theta[l_max]
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
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================
void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
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
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

