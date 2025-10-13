#include <TMB.hpp>
#include <fftw3.h>
#include <unsupported/Eigen/FFT>
#include "../inst/include/predation.hpp"

using namespace atomic;

template<class Type>
Type objective_function<Type>::operator() ()
{ 
  Type nll = 0.0;
  return nll;
}
