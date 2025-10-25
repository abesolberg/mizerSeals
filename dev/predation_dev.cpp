#ifndef EIGEN_FFTW_DEFAULT
#define EIGEN_FFTW_DEFAULT
#endif

#include <TMB.hpp>
#include <complex>
#include <Eigen/Dense>
#include <fftw3.h> // location of FFTW files
#include <unsupported/Eigen/FFT>
#include <fft.hpp>
#include <vector>
#include <complex>
#include <cassert>
#include <cppad/cppad.hpp>
#include <iostream>
#include "../inst/include/predation.h"
#include "../inst/include/helpers.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{ 
  // =========================== DATA ===============================
  
  // Initial Values
  DATA_MATRIX(n) ;
  DATA_VECTOR(n_pp) ;
  
  // Bins
  DATA_VECTOR(w_full) ;
  DATA_VECTOR(dw_full) ;
  DATA_VECTOR(w) ;
  DATA_VECTOR(dw) ;
  DATA_MATRIX(ft_mask) ;
  
  // Pred Kernel
  DATA_MATRIX(ft_pred_kernel_p_real) ;
  DATA_MATRIX(ft_pred_kernel_p_imag) ;
  DATA_MATRIX(ft_pred_kernel_real) ;
  DATA_MATRIX(ft_pred_kernel_imag) ;
  
  // Predation Data Values
  DATA_MATRIX(search_vol) ;
  DATA_MATRIX(intake_max) ;
  
  // Interactions
  DATA_VECTOR(species_interaction_resource) ;
  DATA_MATRIX(species_interaction) ;
  
  // FFT Dims
  DATA_INTEGER(nRows);
  DATA_INTEGER(nCols);
  
  // =========================== Parameters ===============================
  
  PARAMETER(A);
  
  // =========================== Set Values ===============================
  
  std::cout << "*** Setting FFT Dims "  << std::endl;
  set_fft_dims((int)nRows,(int)nCols);
  
  // =========================== Set Bins ===============================
  
  // Calculate Encounter
  std::cout << "*** ENCOUNTER "  << std::endl;
  matrix<Type> ENCOUNTER = getEncounter(species_interaction_resource, n_pp, n, species_interaction, w_full, w , dw_full , ft_pred_kernel_real , ft_pred_kernel_imag, search_vol) ;
  
  // Calculate Feeding Level
  std::cout << "*** FEEDING LEVEL "  << std::endl;
  matrix<Type> FEEDING_LEVEL = getFeedingLevel(intake_max , ENCOUNTER) ;
  
  // Calculate Pred Rate
  std::cout << "*** PRED RATE "  << std::endl;
  matrix<Type> PRED_RATE = getPredRate(n , w_full , w , FEEDING_LEVEL , search_vol , dw , ft_pred_kernel_p_real , ft_pred_kernel_p_imag , ft_mask) ;
  
  // Calculate Pred Mort
  std::cout << "*** PRED MORT "  << std::endl;
  matrix<Type> PRED_MORT = getPredMort(w_full, w, species_interaction, PRED_RATE) ;
  
  // Calculate Resource Mort
  std::cout << "*** RESOURCE MORT "  << std::endl;
  vector <Type> RESOURCE_MORT = getResourceMort(species_interaction_resource , PRED_RATE) ;
  
  // =========================== Calculating Likelihood ===============================
  
  Type nll = 0.0;
  
  // Reports
  REPORT(ENCOUNTER) ;
  REPORT(FEEDING_LEVEL) ;
  REPORT(PRED_RATE) ;
  REPORT(PRED_MORT) ;
  REPORT(RESOURCE_MORT) ;
  
  return nll;
}
