#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  // Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.67;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  // Recombination parameters
  // double Yp          = 0.245;
  double Yp          = 0.0;


  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  bool m1_output = false;
  bool supernova = false;
  bool m2_output = false;
  bool m3_output = false;
  bool m4_output = false;


  //=========================================================================
  // Module I
  //=========================================================================

  // m1_output=true;
  // supernova=true;
  if(m1_output){
    // Set up and solve the background
    std::string M1_OUTPATH = "milestone1/data/";
    
    BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
    cosmo.solve();
    cosmo.info();
    
    if(supernova){
      // Output background evolution quantities
      Utils::StartTiming("Supernova");
      Vector best_fit_params = mcmc_fit_to_supernova_data(
        M1_OUTPATH + "supernovadata.txt", 
        M1_OUTPATH + "supernovafitnew.txt");
      Utils::EndTiming("Supernova"); 
      
      //=====================
      // Run simulation with parameters 
      // from the best supernova fit 
      //=====================
      double h_est          = best_fit_params[0];//0.70189;
      double OmegaM_est     = best_fit_params[1];//0.25932;
      double OmegaK_est     = best_fit_params[2];//0.0673887;
      double OmegaCDM_est   = OmegaM_est - OmegaB; 

      BackgroundCosmology bestSNfit(h_est, OmegaB, OmegaCDM_est, OmegaK_est, 0, TCMB);

      bestSNfit.solve();
      bestSNfit.output_dL(
        M1_OUTPATH + "bestfit_dL.txt", 
        -2.0, 0.0
      );
      cosmo.output_dL(
        M1_OUTPATH + "planck_dL.txt",
        -2.0, 0.0
      );

    }

    else{
      // cosmo.output(M1_OUTPATH + "NEWcosmology_compare.txt");          // 
      cosmo.output(M1_OUTPATH + "cosmology_new.txt");          // Consistency checks and analysis  
      // cosmo.output(M1_OUTPATH + "NEWcosmology_times.txt");    // High resolution for important times  
    }



  }

  //=========================================================================
  // Module II
  //=========================================================================
  if(m2_output){
    // Set up and solve the background
    BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
    cosmo.solve();

    // Solve the recombination history
    RecombinationHistory rec(&cosmo, Yp);
    rec.solve();
    rec.info();
  
  
    // Output recombination quantities
    rec.output("milestone2/data/recombination.txt");
    rec.output("milestone2/data/recombination_saha.txt");
    rec.output("milestone2/data/recombination_split.txt");
  
    //===== Find decoupling times ======
    rec.output_important_times("milestone2/data/rec_times.txt");
    rec.output_important_times("milestone2/data/rec_times_saha.txt");
  }

  //=========================================================================
  // Module III
  //=========================================================================
  if(m3_output){
    // Set up and solve the background
    // Set Neff = 0 in this milestone to ignore neutrinos.
    BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, 0.0, TCMB);
    cosmo.solve();

    // Solve the recombination history
    RecombinationHistory rec(&cosmo, Yp);
    rec.solve();
  
    // Solve the perturbations
    Perturbations pert(&cosmo, &rec);
    pert.solve();
    pert.info();
    
    // Output perturbation quantities
    double kvalue = 0.001 / Constants.Mpc;
    pert.output(kvalue, "milestone3/data/perturbations_k0.001.txt");

    kvalue = 0.03 / Constants.Mpc;
    pert.output(kvalue, "milestone3/data/perturbations_k0.03.txt");

    kvalue = 0.3 / Constants.Mpc;
    pert.output(kvalue, "milestone3/data/perturbations_k0.3.txt");
  }


  m4_output=true;
  //=========================================================================
  // Module IV
  //=========================================================================
  if(m4_output){
    bool source = true;
    // Set up and solve the background
    // Set Neff = 0 in this milestone to ignore neutrinos.
    BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, 0.0, TCMB);
    cosmo.solve();

    // Solve the recombination history
    RecombinationHistory rec(&cosmo, Yp);
    rec.solve();
  
    // Solve the perturbations
    Perturbations pert(&cosmo, &rec);
    pert.solve(source);

    PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
    power.solve(false);

    std::string M4_DATA_PATH = "milestone4/data/";

    // power.output(M4_DATA_PATH + "cells.txt");
    // power.outputPS(M4_DATA_PATH + "matterPS.txt", 1000);
    power.outputThetas(M4_DATA_PATH + "thetas.txt", 2000);
  
  // Remove when module is completed
  }
  return 0;

  // Utils::EndTiming("Everything");
}
