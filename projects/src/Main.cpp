#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

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
  double Yp          = 0.245;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  // cosmo.output("milestone1/data/cosmology.txt");          // Consistency checks and analysis  
  // cosmo.output("milestone1/data/cosmology_dL.txt");       // Comparing with supernova data 
  // cosmo.output("milestone1/data/cosmology_times.txt");    // High resolution for important times  

  //=====================
  // Run simulation with parameters 
  // from the best supernova fit 
  //=====================
  // double h_est          = 0.70189;
  // double OmegaM_est     = 0.25932;
  // double OmegaK_est     = 0.0673887;
  // double OmegaCDM_est   = OmegaM_est - OmegaB;
  // BackgroundCosmology bestSNfit(h_est, OmegaB, OmegaCDM_est, OmegaK_est, Neff, TCMB);
  // bestSNfit.solve();
  // bestSNfit.info();
  // bestSNfit.output("milestone1/data/bestfit_cosmology_dL.txt");

  // Utils::StartTiming("Supernova");
  // mcmc_fit_to_supernova_data("milestone1/data/supernovadata.txt", "milestone1/data/supernovafit.txt");
  // Utils::EndTiming("Supernova"); 

  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");
  
  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.01 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.01.txt");
  
  // Remove when module is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
