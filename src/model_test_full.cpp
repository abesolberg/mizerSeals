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


// params values removed 
// -- intake_max
// -- species_params_alpha
// -- metab
// -- psi
// -- w_full
// -- w
// -- selectivity
// -- catchability
// -- mu_b
// -- dw
// -- species_params_erepro
// -- w_min_idx
// -- interaction
// -- species_params_rmax
// -- species_params_interaction_resource
// -- species_params_interaction
// -- rr_pp
// -- cc_pp
// -- dt
// -- dw_full
// -- ft_pred_kernel_e
// -- search_vol
// -- ft_pred_kernel_p
// -- ft_mask

// needs initial conditions
// -- n_pp
// -- n

// Effort required to run dynamics


// Global storage for dimensions
static int g_nRows = 0;  // Default value
static int g_nCols = 0;  // Default value

// Function to set nRows and nCols before calling the atomic function
void set_fft_dims(int nR, int nC) {
  g_nRows = nR;
  g_nCols = nC;
}


//// useful functions
///  can't find this function so making it
template<class Type>
vector<Type> vec_MS(matrix<Type> M_in){
  int n_rows = M_in.rows();
  vector<Type> ret(n_rows * M_in.cols());
  for (int i = 0; i < M_in.cols(); i++) {
    ret.segment(i * n_rows,n_rows) = M_in.col(i);
  }
  return(ret);
}

/// exponential of a matrix element by element
template<class Type>
matrix<Type> exp_mat_ele(matrix<Type> M_in){
  matrix<Type> ret = M_in;
  for (int i = 0; i < M_in.rows(); i++) {
    for (int j = 0; j < M_in.cols(); j++) {
      ret(i,j) = exp(M_in(i,j));
    }
  }
  return(ret);
}

// make everything in an array an exponential (only works for 2D array)
template<class Type>
array<Type> exp_arr_ele(array<Type> M_in){
  array<Type> ret = M_in;
  for (int i = 0; i < M_in.rows(); i++) {
    for (int j = 0; j < M_in.cols(); j++) {
      ret(i,j) = exp(M_in(i,j));
    }
  }
  return(ret);
}

// Set vector as diagonal of matrix -- would like to remove dummy matrix but can't figure it out
template<class Type>
array<Type> vec_to_diag(array<Type>dummy , vector<Type> vec){
  array<Type> ret = dummy ;
  for (int i = 0; i < ret.rows(); i++) {
    ret(i,i) = vec(i);
  }
  return(ret);
}


//// functions involving FFT
template <typename T>
inline double get_base(const T &x) { 
  return x; 
}
template <typename T>
inline double get_base(const CppAD::AD<T>& x) { 
  return get_base( CppAD::Value(x) );
}


template<class Type>
matrix<Type> computePrey(vector<Type> species_params_interaction_resource, 
                         vector<Type> n_pp, 
                         matrix<Type> n, 
                         matrix<Type> species_params_interaction, 
                         vector<Type> w_full, 
                         vector<Type> w,
                         vector<Type> dw_full) {
  matrix<Type> ret(species_params_interaction_resource.size(), n_pp.size());
  for (int i = 0; i < ret.rows(); ++i) {
    for (int j = 0; j < ret.cols(); ++j) {
      ret(i,j) = species_params_interaction_resource(i) * n_pp(j);
    }
  }
  ret.block(0, w_full.size() - w.size(), ret.rows(), w.size()) += species_params_interaction * n;
  vector<Type> tmp = w_full * dw_full;
  for (int i = 0; i < ret.rows(); ++i) {
    ret.row(i) = vector<Type>(ret.row(i)) * tmp;
  }
  return ret;
}


template <class Type>
Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> 
convertToRowMajor(const matrix<Type>& M) {
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixRM;
  MatrixRM M_rm(M.rows(), M.cols());
  for (int i = 0; i < M.rows(); ++i) {
    for (int j = 0; j < M.cols(); ++j) {
      M_rm(i, j) = M(i, j);
    }
  }
  return M_rm;
}


template <class Type>
void buildComplexKernel(
    const matrix<Type>& RealPart, 
    const matrix<Type>& ImagPart,
    matrix<Type>& kernel_real, 
    matrix<Type>& kernel_imag) 
{
  // Resize before assigning values
  kernel_real.resize(RealPart.rows(), RealPart.cols());
  kernel_imag.resize(ImagPart.rows(), ImagPart.cols());
  
  // Assign values manually (avoiding eval())
  for (int i = 0; i < RealPart.rows(); ++i) {
    for (int j = 0; j < RealPart.cols(); ++j) {
      kernel_real(i, j) = RealPart(i, j);
      kernel_imag(i, j) = ImagPart(i, j);
    }
  }
}


template <class Type>
vector<std::complex<Type> > prepareFFTInput(const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& M) {
  int nRows = M.rows();
  int nCols = M.cols();
  vector<std::complex<Type> >fft_in(nRows * nCols);
  
  for (int i = 0; i < nRows; i++){
    for (int j = 0; j < nCols; j++){
      int idx = (i * nCols + j);
      fft_in[idx] = std::complex<Type>(M(i,j), Type(0));
    }
  }
  return fft_in;
}


namespace atomic {

// Row-wise FFT: process each row separately
// (Assuming this function is declared as a template so that 'inverse' is a compile-time constant.)
template<bool inverse = false>
void batched_fft_work(const CppAD::vector<double>& x,
                      CppAD::vector<double>& y)
{
  // Read dimensions from the first two entries.
  int nRows = g_nRows;
  int nCols = g_nCols;
  
  y.resize(x.size());
  
  typedef std::complex<double> C;
  Eigen::FFT<double> f;
  f.SetFlag(f.Unscaled);
  
  // Loop over each row.
  for (int r = 0;  r < nRows; r++) {
    // Compute the offset (in doubles) for row r.
    int offset = 2 * (r) * nCols;
    // does this casting do what i want it do do??
    
    C* X = (C*)(x.data() + offset);
    C* Y = (C*)(y.data() + offset);
    // Specify the number of complex numbers in this row.
    int ncplx = nCols;
    if (!inverse)
      f.fwd(Y, X, ncplx);
    else
      f.inv(Y, X, ncplx);
  }
  
}


// Declare the atomic functions.
TMB_ATOMIC_VECTOR_FUNCTION_DECLARE(batched_fft)
  TMB_ATOMIC_VECTOR_FUNCTION_DECLARE(batched_ifft)
  
  // Define the forward atomic function.
  TMB_ATOMIC_VECTOR_FUNCTION_DEFINE(
    batched_fft,
    tx.size(),
    batched_fft_work<0>(tx, ty),  // forward pass (FFT)
    px = batched_ifft(py)             // reverse pass calls iFFT
  )
  
  // Define the inverse atomic function.
  TMB_ATOMIC_VECTOR_FUNCTION_DEFINE(
    batched_ifft,
    tx.size(),
    batched_fft_work<1>(tx, ty),   // forward pass (iFFT)
    px = batched_fft(py)              // reverse pass calls FFT
  )
  
  // Helper: convert a standard vector to a CppAD vector, call the atomic operator, and convert back.
  template<class Type>
  vector<std::complex<Type> > batched_fft_b(vector<std::complex<Type> > x,
                                            bool inverse = false) {
    CppAD::vector<Type> x_ad(2 * x.size());
    Type* px = x_ad.data();
    Type* pxc = (Type*) x.data();
    for (size_t i=0; i<x_ad.size(); i++) {
      px[i] = pxc[i];
    }
    
    CppAD::vector<Type> y_ad = inverse ? atomic::batched_ifft(x_ad) : atomic::batched_fft(x_ad);
    px = (Type*) y_ad.data();
    for (size_t i=0; i<x_ad.size(); i++) {
      pxc[i] = px[i];
    }
    return x;
  }
  
}

template <class Type>
vector<std::complex<Type>> scaleFFTOutput(const vector<std::complex<Type>> &fft_data, int nCols) {
  vector<std::complex<Type>> scaled(fft_data.size());  // Correctly initialize output vector
  
  Type scale = static_cast<Type>(nCols);  // Ensure correct type
  
  for (int i = 0; i < fft_data.size(); i++) {
    scaled[i] = fft_data[i] / scale;  // Scale both real and imaginary parts
  }
  return scaled;
}


template <class Type>
vector<std::complex<Type>> multiplyByKernel(
    vector<std::complex<Type>>& fft_data,  // Pass-by-reference for efficiency
    const matrix<Type>& kernel_real, 
    const matrix<Type>& kernel_imag, 
    int nRows, int nCols) 
{
  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      std::complex<Type> kernel_val(kernel_real(i, j), kernel_imag(i, j));
      
      // Apply kernel multiplication to the correct index
      int idx = i * nCols + j;  // Assuming row-major order
      fft_data[idx] *= kernel_val;
    }
  }
  
  return fft_data;  // Missing return fixed
}



template <class Type>
Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
convertFFTOutputToMatrix(const std::vector<Type>& fft_data, int nRows, int nCols) {
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixRM;
  MatrixRM M(nRows, nCols);
  for (int i = 0; i < nRows; i++){
    for (int j = 0; j < nCols; j++){
      int idx = (i * nCols + j);
      M(i,j) = static_cast<Type>(fft_data[idx]);
    }
  }
  return M;
}

template <class Type>
matrix<Type> computeQMatrix(
    matrix<Type> n, 
    vector<Type> w_full, 
    vector<Type> w, 
    matrix<Type> feeding_level, 
    matrix<Type> search_vol, 
    vector<Type> dw
) 
{
  matrix<Type> Q_mat(n.rows(), w_full.size());
  Q_mat.setZero();
  
  matrix<Type> tmp = (Type(1.0) - feeding_level.array()) * n.array() * search_vol.array();
  tmp = tmp.array().rowwise() * dw.transpose().array();
  
  // Place tmp into the right-hand block of Q_mat
  Q_mat.block(0, w_full.size() - w.size(), Q_mat.rows(), w.size()) = tmp;
  return Q_mat;
}

template<class T>
vector<std::complex<T> > cplx(vector<T> x) {
  vector<std::complex<T> > xc(x.size());
  for (int i = 0; i < x.size(); i++) {
    xc[i] = std::complex<T>(x[i], T(0));
  }
  return xc;
}


template <class Type>
matrix<Type> get_pred_rate(
    matrix<Type> n, 
    vector<Type> w_full, 
    vector<Type> w, 
    matrix<Type> feeding_level, 
    matrix<Type> search_vol, 
    vector<Type> dw, 
    matrix<Type> ft_pred_kernel_p_real,
    matrix<Type> ft_pred_kernel_p_imag, 
    matrix<Type> ft_mask
) {
  
  matrix<Type> Q_mat = computeQMatrix(n, w_full, w, feeding_level, search_vol, dw);
  // Convert to row-major
  auto Q_rm = convertToRowMajor(Q_mat);
  int nRows = Q_rm.rows();
  int nCols = Q_rm.cols();
  // Build the complex kernel for pred rate.
  matrix<Type> kernel_real, kernel_imag;
  buildComplexKernel(ft_pred_kernel_p_real, ft_pred_kernel_p_imag, kernel_real, kernel_imag);
  
  // Prepare FFTW input
  vector<std::complex<Type> > fft_in = prepareFFTInput(Q_rm);
  
  auto fft_forward = atomic::batched_fft_b(fft_in);
  
  auto fft_forward_scaled = scaleFFTOutput(fft_forward, nCols);
  
  auto mult = multiplyByKernel(fft_forward_scaled, kernel_real, kernel_imag, nRows, nCols);
  
  vector<Type> fft_inverse = atomic::batched_fft_b(mult, true).real();
  
  auto pred_matrix = convertFFTOutputToMatrix<Type>(fft_inverse, nRows, nCols);
  
  auto ft_mask_rm = convertToRowMajor(ft_mask);
  
  auto ret_matrix = pred_matrix.array() * ft_mask_rm.array();
  
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixRM;
  
  // Force negative values to zero.
  MatrixRM ret_matrix_mutable = ret_matrix;
  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      if (ret_matrix_mutable(i, j) < Type(0.0)) {
        ret_matrix_mutable(i, j) = Type(0.0);
      }
    }
  }
  //can i just return directly?
  matrix<Type> ret = ret_matrix_mutable;
  return ret;
}


template <class Type>
matrix<Type> getEncounter(
    vector<Type> species_params_interaction_resource, 
    vector<Type> n_pp, 
    matrix<Type> n, 
    matrix<Type> species_params_interaction, 
    vector<Type> w_full, 
    vector<Type> w, 
    vector<Type> dw_full, 
    matrix<Type> ft_pred_kernel_real,
    matrix<Type> ft_pred_kernel_imag, 
    matrix<Type> search_vol
) {
  // Compute prey matrix.
  matrix<Type> prey = computePrey(species_params_interaction_resource, n_pp, n, species_params_interaction, w_full, w, dw_full);
  
  
  // Convert to row-major.
  auto prey_rm = convertToRowMajor(prey);
  
  int nRows = prey_rm.rows();
  int nCols = prey_rm.cols();
  
  matrix<Type> kernel_real, kernel_imag;
  buildComplexKernel(ft_pred_kernel_real, ft_pred_kernel_imag, kernel_real, kernel_imag);
  
  vector<std::complex<Type> > fft_in = prepareFFTInput(prey_rm);
  
  auto fft_forward = atomic::batched_fft_b(fft_in);
  
  auto fft_forward_scaled = scaleFFTOutput(fft_forward, nCols);
  
  auto mult = multiplyByKernel(fft_forward_scaled, kernel_real, kernel_imag, nRows, nCols);
  
  vector<Type> fft_inverse = atomic::batched_fft_b(mult, true).real();
  
  auto encounter_matrix = convertFFTOutputToMatrix<Type>(fft_inverse, nRows, nCols);
  
  
  // Extract desired block: columns (w_full.size()-w.size()) to end.
  int start_col = static_cast<int>(w_full.size()) - static_cast<int>(w.size());
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixRM;
  MatrixRM ret_block = encounter_matrix.block(0, start_col, nRows, w.size());
  
  
  // GETDIET OUTPUT GOES TO HERE 
  
  //THIS IS WHAT THE GETDIETENCOUNTER OUTPUT WIL GET TO.
  //inter <- cbind(params@interaction, params@species_params$interaction_resource)
  //diet[, , 1:(no_sp + 1)] <- sweep(sweep(diet[, , 1:(no_sp + 
  //1), drop = FALSE], c(1, 3), inter, "*"), c(1, 2), params@search_vol, 
  //"*")
  //make the getdietencounter function that runs every 10 timesteps.
  //the getdiet function will be the same as the entire getdiet function from mizer 
  //till the line above.
  
  //GET THE DIET FUNCTION WORKING SO THAT IT OUTPUTS DIET ALSO
  //
  
  
  // Multiply elementwise by search_vol.
  ret_block = ret_block.array() * convertToRowMajor(search_vol).array();
  
  matrix<Type> ret = ret_block;
  return ret;
}



/// other rates functions

template<class Type>
matrix <Type> FeedingLevel(matrix<Type> intake_max, matrix<Type> encounter){
  matrix<Type> ret(encounter.rows(),encounter.cols()) ;
  for (int i=0; i<ret.rows(); ++i){
    for (int j=0; j<ret.cols(); ++j){
      ret(i,j) = encounter(i,j) / (encounter(i,j) + intake_max(i,j));
    }
  }
  return ret;
}

// Metabolism
template<class Type> 
matrix <Type> getMetab(vector<Type> w , vector<Type> k , vector<Type> ks , vector<Type> p){
  matrix<Type> ret(k.size() , w.size()) ;
  for (int i=0; i<ret.rows(); ++i){
    for (int j=0; j<ret.cols(); ++j){
      ret(i,j) = pow(w(j),p(i))*ks(i) + (k(i)*w(j));
    }
  }
  return ret;
}

template<class Type>
matrix <Type> EReproAndGrowth(matrix<Type> feeding_level, matrix<Type> encounter, vector<Type> species_params_alpha, matrix<Type> metab){
  matrix<Type> ret(feeding_level.rows(),feeding_level.cols()) ;
  for (int i=0; i<ret.rows(); ++i){
    for (int j=0; j<ret.cols(); ++j){
      ret(i,j) = (Type(1.0) - feeding_level(i,j)) * encounter(i,j) * species_params_alpha(i) - metab(i,j);
    }
  }
  return ret;
}

template<class Type>
matrix <Type> ERepro(matrix<Type> psi, matrix<Type> erepog){
  matrix<Type> ret(psi.rows(),psi.cols()) ;
  for (int i=0; i<ret.rows(); ++i){
    for (int j=0; j<ret.cols(); ++j){
      if (erepog(i,j) < Type(0.0)){ // is this required?
        erepog(i,j) = Type(0.0);
      }
      ret(i,j) = erepog(i,j) * psi(i,j);
    }
  }
  return ret;
}

template<class Type>
matrix <Type> EGrowth(matrix<Type> erepog, matrix<Type> e_repro){
  matrix<Type> ret(erepog.rows(),erepog.cols()) ;
  for (int i=0; i<ret.rows(); ++i){
    for (int j=0; j<ret.cols(); ++j){
      if (erepog(i,j) < Type(0.0)){ // is this required?
        erepog(i,j) = Type(0.0);
      }
    }
  }
  ret = erepog - e_repro;
  return ret;
}

template<class Type>
matrix <Type> PredMort(vector<Type> w_full, vector<Type> w, matrix<Type> interaction, matrix<Type> pred_rate){
  matrix<Type> ret(pred_rate.rows(),w.size()) ;
  matrix <Type> tmp= pred_rate.block(0,w_full.size()-w.size(),pred_rate.rows(),w.size());
  ret = interaction.transpose() * tmp; // ask Gustav about drop==FALSE
  return ret;
}

// this one is within computational precision -- may need experiment and change
template<class Type>
vector <Type> ResourceMort(vector<Type> species_params_interaction_resource, matrix<Type> pred_rate){
  vector<Type> ret(pred_rate.cols()) ;
  matrix<Type> tmp = pred_rate.transpose();
  ret = tmp * species_params_interaction_resource; 
  return ret;
}


template<class Type>
array<Type> mizerFMortGear(array<Type> selectivity, array<Type> catchability, vector<Type> effort){
  array<Type> ret = selectivity;
  array<Type> tmp = catchability.transpose();
  for (int i=0; i<catchability.rows(); ++i){
    tmp.col(i) = tmp.col(i) * effort(i);
  }
  int tmp1 = ret.cols();
  tmp = tmp.transpose();
  
  for (int i=0; i < tmp1; ++i){
    ret.col(i) = ret.col(i) * tmp;
  }
  
  
  return ret;
}

template<class Type>
array<Type> FMort(array<Type> selectivity, array<Type> catchability, vector<Type> effort){
  array<Type> ret(catchability.cols(),selectivity.cols());
  array<Type> tmp = mizerFMortGear(selectivity, catchability, effort);
  
  for (int i=0; i<ret.rows(); ++i){
    for (int j=0; j<ret.cols(); ++j){
      ret(i,j) = sum(tmp.col(j).col(i));
      //ret(i,j) = tmp(i,1,j);
    }
  }
  return ret;
}

// Z Mort Size
template<class Type>
array<Type> mizerZMortSize(array<Type> Zselectivity, array<Type> Zcatchability, vector<Type> Zscale){
  array<Type> ret = Zselectivity;
  array<Type> tmp = Zcatchability.transpose();
  for (int i=0; i<Zcatchability.rows(); ++i){ //it may be worth only scaling z > 0 
    tmp.col(i) = tmp.col(i) * Zscale(i);
  }
  int tmp1 = ret.cols();
  tmp = tmp.transpose();
  
  for (int i=0; i < tmp1; ++i){
    ret.col(i) = ret.col(i) * tmp;
  }
  
  
  return ret;
}

// Z Mort
template<class Type>
array<Type> ZMort(array<Type> Zselectivity, array<Type> Zcatchability, vector<Type> Zscale){
  array<Type> ret(Zcatchability.cols(),Zselectivity.cols());
  array<Type> tmp = mizerZMortSize(Zselectivity, Zcatchability, Zscale);
  
  for (int i=0; i<ret.rows(); ++i){
    for (int j=0; j<ret.cols(); ++j){
      ret(i,j) = sum(tmp.col(j).col(i));
      //ret(i,j) = tmp(i,1,j);
    }
  }
  return ret;
}

template<class Type>
matrix<Type> Mort(matrix<Type> pred_mort, matrix<Type> mu_b, array<Type> f_mort , array<Type> z_mort){
  matrix<Type> ret = pred_mort + mu_b + f_mort.matrix() + z_mort.matrix();
  // no other mortality at the moment
  return ret;
}

// template<class Type>
// matrix<Type> Mort(matrix<Type> pred_mort, matrix<Type> mu_b, array<Type> f_mort){
//   matrix<Type> ret = pred_mort + mu_b + f_mort.matrix();
//   // no other mortality at the moment
//   return ret;
// }

// template<class Type>
// vector<Type> RDI(matrix<Type> e_repro, matrix<Type> n, vector<Type> dw, vector<Type> species_params_erepro, vector<Type> w, vector<int> w_min_idx){
//   vector <Type> ret(e_repro.rows());
//   matrix <Type> tmp1 = e_repro.array() * n.array();
//   vector <Type> tmp2 = tmp1 * dw;
//   for (int i=0; i<ret.size(); ++i){
//     ret(i) = Type(0.5) * (tmp2(i) * species_params_erepro(i)) / w(w_min_idx(i)-int(1)); // less 1 because its cpp and not R. 
//   }
//   return ret;
// }

template<class Type>
vector<Type> RDI(matrix<Type> e_repro, matrix<Type> n, vector<Type> dw, vector<Type> species_params_erepro, vector<Type> w, vector<int> w_min_idx , vector<Type> female_ratio){
  vector <Type> ret(e_repro.rows());
  matrix <Type> tmp1 = e_repro.array() * n.array();
  vector <Type> tmp2 = tmp1 * dw;
  for (int i=0; i<ret.size(); ++i){
    ret(i) = female_ratio(i) * (tmp2(i) * species_params_erepro(i)) / w(w_min_idx(i)-int(1)); // less 1 because its cpp and not R. Added sex ratio here. Should be .5, 1, or 0
  }
  return ret;
}

template<class Type>
matrix<Type> RDD(vector<Type> rdi, vector<Type> species_params_rmax){
  vector<Type> ret(rdi.size());
  for (int i=0; i<ret.size(); ++i){
    ret(i) = rdi(i) / (Type(1.0) + (rdi(i) / species_params_rmax(i))); 
  }
  return ret;
}

template<class Type>
matrix<Type> RDD2(vector<Type> rdi, vector<Type> species_params_rmax , vector<Type> RScale, vector<int> repro_idx , vector<Type> sex_split){
  vector<Type> ret(rdi.size());
  for (int i=0; i<ret.size(); ++i){
    ret(i) = (rdi(i) / (Type(1.0) + (rdi(i) / species_params_rmax(i))))*RScale(i); 
  }
  for (int i=0; i<ret.size(); ++i){
    ret(i) = ret(repro_idx(i)-int(1))/sex_split(i); //repro idx should be the index for female RDI, sex split = 1 unless repro index is not unique, then = 2
  }
  
  return ret;
}

template<class Type>
array <Type> set_allo_mort(vector<Type>w , vector<Type> z0pre , vector<Type> z0exp) {
  array<Type> ret(z0pre.size() , z0pre.size() , w.size()) ;
  ret.setZero() ;
  for(int i=0 ; i < ret.rows() ; ++i) {
    for(int j=0 ; j < ret.cols() ; ++j) {
      ret(i,i,j) = z0pre(i)*pow(w(j),-z0exp(i)) ;
    }
  }
  return ret ;
}

template<class Type>
vector <Type> resource_semichemostat(vector<Type> rr_pp, vector<Type> resource_mort, vector<Type> cc_pp, vector<Type> n_pp, Type dt){
  //mur <- resource_rate + rates$resource_mort
  vector<Type> mur = rr_pp + resource_mort;
  //n_steady <- resource_rate * resource_capacity/mur
  vector<Type> n_steady = rr_pp.array() * cc_pp.array()/mur.array();
  //n_pp_new <- n_steady + (n_pp - n_steady) * exp(-mur * dt)
  vector<Type> n_pp_new = n_steady + (n_pp - n_steady) * exp(-mur * dt);
  //sel <- !is.finite(n_pp_new)
  //n_pp_new[sel] <- n_pp[sel]
  //n_pp_new
  return n_pp_new;
}

template<class Type>
matrix<Type> getA(matrix <Type> e_growth, vector<Type> dw, Type dt){
  matrix<Type> ret(e_growth.rows(),e_growth.cols());
  ret.setZero();
  matrix<Type> tmp = e_growth.block(0,0,ret.rows(),ret.cols()-1) * dt;
  vector<Type> tele = dw.segment(1,dw.size()-1);
  for (int i=0; i<tmp.rows(); ++i){
    for (int j=0; j<tele.size(); ++j){
      tmp(i,j) = -tmp(i,j) / tele(j); 
    }
  }
  ret.block(0,1,ret.rows(),ret.cols()-1) = tmp;
  return ret;
}

template<class Type>
matrix<Type> getB(matrix <Type> e_growth, vector<Type> dw, Type dt,matrix <Type> mort){
  matrix<Type> ret(e_growth.rows(),e_growth.cols());
  for (int i=0; i<ret.rows(); ++i){
    for (int j=0; j<dw.size(); ++j){
      ret(i,j) = 1 + e_growth(i,j) * dt / dw(j) + mort(i,j) * dt; 
    }
  }
  return ret;
}

template<class Type>
matrix<Type> getS(matrix<Type> n){
  matrix<Type> ret(n.rows(),n.cols());
  ret.setZero();
  ret.block(0,1,ret.rows(),ret.cols() - 1) = n.block(0,1,ret.rows(),ret.cols() - 1);
  return ret;
}

template <class Type>
matrix<Type> inner_project_loop(matrix<Type> n, matrix<Type> A, matrix<Type> B, matrix<Type> S, vector<int> w_min_idx){
  int no_sp = n.rows();
  int no_w = n.cols();
  for (int i = 0; i < no_sp; i++) {
    for (int j = w_min_idx[i]; j < no_w; j++) {
      n(i,j) = (S(i,j) - A(i,j)*n(i,j-1)) / B(i,j);
    }
  }
  return n;
}


template<class Type>
matrix<Type> get_N_after(vector<Type> vec_in, int n_rows, int n_cols){
  matrix<Type> ret(n_rows, n_cols);
  for (int i = 0; i < n_cols; i++) {
    ret.col(i) = vec_in.segment(i * n_rows,n_rows);
  }
  return(ret);
}

template<class Type>
vector<Type> get_Npp_after(vector<Type> vec_in, int n_rows, int n_cols){
  vector<Type> ret(vec_in.size() - n_rows * n_cols);
  ret = vec_in.segment(n_cols * n_rows,ret.size());
  return(ret);
}


template<class Type>
matrix<Type> running_model(vector<Type> species_params_interaction_resource, vector<Type> n_pp, matrix<Type> n, matrix <Type> species_params_interaction, vector<Type> w_full, vector<Type> w,vector<Type> dw_full,
                           matrix<Type> ft_pred_kernel_real, matrix<Type> ft_pred_kernel_imag, matrix<Type> search_vol,matrix<Type> intake_max,vector<Type> species_params_alpha, //matrix<Type> metab,
                           matrix<Type> psi, vector<Type> dw, matrix<Type> ft_pred_kernel_p_real, matrix<Type> ft_pred_kernel_p_imag, matrix<Type> ft_mask, array<Type> f_mort, matrix<Type> mu_b, 
                           vector<Type> species_params_erepro, vector<int> w_min_idx, vector<Type> species_params_rmax, vector<Type> rr_pp, vector<Type> cc_pp, 
                           Type dt , 
                           vector<Type> k , vector<Type> ks , vector<Type> p ,
                           array<Type> z_mort , array<Type> trawl_catch , 
                           vector<int> repro_idx, vector<Type> sex_split, vector<Type> female_ratio ,
                           vector<Type> r_scale){
  
  // in a loop calculate
  matrix <Type> encounter = getEncounter(species_params_interaction_resource, n_pp, n, species_params_interaction, w_full, w, dw_full,ft_pred_kernel_real, ft_pred_kernel_imag,search_vol);
  matrix <Type> feeding_level = FeedingLevel(intake_max, encounter);
  matrix <Type> metab = getMetab(w , k , ks , p) ;
  matrix <Type> erepog = EReproAndGrowth(feeding_level, encounter, species_params_alpha, metab);// is e
  matrix <Type> e_repro = ERepro(psi,erepog);
  matrix <Type> e_growth = EGrowth(erepog,e_repro);
  matrix <Type> pred_rate = get_pred_rate(n, w_full, w, feeding_level, search_vol, dw, ft_pred_kernel_p_real, ft_pred_kernel_p_imag, ft_mask);
  matrix <Type> pred_mort = PredMort(w_full, w, species_params_interaction, pred_rate);
  //array<Type> f_mort = FMort(selectivity, catchability, effort); // this does not depend on the dynamics
  array <Type> selectivity;
  array <Type> catchability;
  array <Type> Tselectivity ;
  array <Type> Tcatchability ;
  array <Type> effort;
  array <Type> trawl_effort ;
  
  matrix<Type> mort = Mort(pred_mort, mu_b, f_mort , z_mort);
  vector<Type> rdi = RDI(e_repro, n , dw , species_params_erepro, w, w_min_idx , female_ratio); // Vector with the rate of egg production for each species
  matrix<Type> rdd = RDD2(rdi,species_params_rmax , r_scale , repro_idx , sex_split); // The flux enetering the smallest size class of each species -- This is recruitment
  vector<Type> resource_mort = ResourceMort(species_params_interaction_resource, pred_rate);
  
  vector<Type> n_pp_new = resource_semichemostat(rr_pp, resource_mort, cc_pp, n_pp, dt);
  // calc new n using A,b,S and the rest of the stuff
  matrix<Type> a = getA(e_growth, dw, dt);
  matrix<Type> b = getB(e_growth, dw, dt,mort);
  matrix<Type> S = getS(n);
  
  matrix<Type> n_new=n;
  
  vector<int> temp_min = w_min_idx - int(1); // sort out the indexing
  /// add recruits
  for (int i=0; i<n.rows(); ++i){
    n_new(i,temp_min(i)) = (n_new(i,temp_min(i)) + rdd(i) * dt / dw(temp_min(i)) ) / b(i,temp_min(i)); 
  }
  
  n_new = inner_project_loop(n_new, a, b, S, w_min_idx);
  
  // end loop
  vector<Type> tmp = vec_MS(n_new);
  vector<Type> ret(n.rows() * n.cols() + n_pp.size());
  ret << tmp,n_pp_new;
  
  return ret;
}

template<class Type>
vector<Type> get_catch(matrix<Type> saved_n, array<Type> f_mort, vector<Type> w_dw, Type dt){
  vector<Type> ret = saved_n.col(0);
  matrix<Type> biomass = saved_n;
  for (int i=0; i<biomass.rows(); ++i){
    biomass.row(i) = vector<Type>(saved_n.row(i)) * w_dw; 
  }
  
  matrix<Type> yield = biomass.array() * f_mort.matrix().array() * dt;
  ret = yield.rowwise().sum();
  
  return(ret);
}

template<class Type>
vector<Type> get_biomass(matrix<Type> saved_n, vector<Type> w_dw, Type dt){
  vector<Type> ret = saved_n.col(0);
  matrix<Type> biomass = saved_n;
  for (int i=0; i<biomass.rows(); ++i){
    biomass.row(i) = vector<Type>(saved_n.row(i)) * w_dw; 
  }
  
  matrix<Type> bm = biomass.array() * dt;
  ret = bm.rowwise().sum();
  
  return(ret);
}

template<class Type>
array<Type> sort_annual_effort(array<Type> selectivity, array<Type> catchability, matrix<Type> effort,int times, Type dt){
  array<Type> ret(catchability.cols(),selectivity.cols(),times);
  int year = 0;
  for (int j=0; j<times; ++j){
    if( (Type(j) * dt) > (Type(year) + 1.0)){
      year += int(1);
    }
    ret.col(j) = FMort(selectivity, catchability, vector<Type>(effort.row(year)));
  }
  return(ret);
}

template<class Type>
array<Type> sort_annual_z(array<Type> Zselectivity, array<Type> Zcatchability, matrix<Type> Zscale,int times, Type dt){
  array<Type> ret(Zcatchability.cols(),Zselectivity.cols(),times);
  int year = 0;
  for (int j=0; j<times; ++j){
    if( (Type(j) * dt) > (Type(year) + 1.0)){
      year += int(1);
    }
    ret.col(j) = ZMort(Zselectivity, Zcatchability, vector<Type>(Zscale.row(year)));
  }
  return(ret);
}

template<class Type>
matrix<Type> sort_annual_r(matrix<Type> Rscale,int times, Type dt){
  matrix<Type> ret(Rscale.cols(),times);
  int year = 0;
  for (int j=0; j<times; ++j){
    if( (Type(j) * dt) > (Type(year) + 1.0)){
      year += int(1);
    }
    ret.col(j) = Rscale.row(year);
  }
  return(ret);
}

template<class Type>
array<Type> get_ann_catches(array<Type> catches, Type dt, int years) {
  array<Type> ret(catches.rows(),years);
  
  int j = 0;
  for (int i=0; i<years; ++i){
    Type t=0.0;
    vector<Type> catch_count(catches.rows());
    catch_count.setZero();
    while(t < 0.999){
      catch_count += catches.col(j);
      t += dt;
      j += 1;
    }
    ret.col(i) = catch_count;
  }
  
  return (ret);
}

template<class Type>
array<Type> get_ann_diet(array<Type> diet , Type dt , int years) {
  array<Type> ret(diet.rows() , diet.rows()+1 , years) ;
  ret.setZero() ;
  int j = 0 ;//# needs to be = 0 in TMB
  for (int i=0 ; i < years ; ++i) {
    Type t = 0.1 ; 
    matrix<Type> a(diet.rows() , diet.rows()+1) ;
    a.setZero() ;
    a = diet.col(j).matrix()*dt ;
    while(t < 0.999) {
      a += (diet.col(j).matrix()*dt) ;
      t += dt;
      j += 1;
    }
    ret.col(i) = a ;
  }
  return (ret) ;
}

template<class Type>
array<Type> get_ann_biomass(array<Type> biomass, Type dt, int years) {
  array<Type> ret(biomass.rows(),years);
  
  int j = 0;
  for (int i=0; i<years; ++i){
    Type t=0.0;
    vector<Type> biomass_count(biomass.rows());
    biomass_count.setZero();
    while(t < 0.999){
      biomass_count += biomass.col(j);
      t += dt;
      j += 1;
    }
    ret.col(i) = biomass_count;
  }
  
  return (ret);
}

// Diet stuff 

template<class Type>
matrix<Type> computePreyDiet(vector<Type> species_params_interaction_resource, 
                             vector<Type> n_pp, 
                             matrix<Type> n, 
                             vector<Type> w_full, 
                             vector<Type> w,
                             vector<Type> dw_full) {
  matrix<Type> ret(species_params_interaction_resource.size()+1, n_pp.size());
  ret.setZero();
  
  for (int i = 0 ; i < ret.cols() ; ++i ) {
    ret(7 , i) = n_pp(i) ;
  }
  
  // resource availability
  
  ret.block(0, w_full.size() - w.size(), 7 , w.size()) = n;
  vector<Type> tmp = w_full * dw_full;
  for (int i = 0; i < ret.rows(); ++i) {
    ret.row(i) = vector<Type>(ret.row(i)) * tmp;
  }
  return ret;
  
}

template <class Type>
array <Type> getDiet(
    vector<Type> species_params_interaction_resource, 
    vector<Type> n_pp, 
    matrix<Type> n, 
    matrix<Type> species_params_interaction, 
    vector<Type> w_full, 
    vector<Type> w, 
    vector<Type> dw_full, 
    matrix<Type> ft_pred_kernel_real,
    matrix<Type> ft_pred_kernel_imag, 
    matrix<Type> search_vol ,
    matrix<Type> feeding_level
) {
  // Compute prey matrix.
  
  matrix<Type> inter(species_params_interaction.rows() , species_params_interaction.rows()+1);
  inter.setZero();
  
  for (int i=0; i<species_params_interaction.rows(); ++i){
    for (int j=0; j<species_params_interaction.cols(); ++j){
      inter(i,j) = species_params_interaction(i,j);
    }
  }
  
  for(int j=0 ; j < species_params_interaction_resource.size() ; ++j) {
    inter(j , 7) = species_params_interaction_resource(j) ; // I would prefer not to hard-code this in, but I can't get it otherwise
  }
  
  matrix<Type> prey = computePreyDiet(species_params_interaction_resource, n_pp, n, w_full, w, dw_full);
  
  array<Type> diet(species_params_interaction.rows() , w.size() , inter.cols()) ;
  diet.setZero() ;
  
  for (int i = 0 ; i < diet.cols() ; ++i ) {
    
    matrix<Type> tmp_prey(inter.rows() , prey.cols()) ;
    
    for (int j = 0 ; j < inter.rows() ; j++) {
      tmp_prey.row(j) = inter(j , i) * prey.row(i) ;
    }
    
    // Convert to row-major.
    auto prey_rm = convertToRowMajor(tmp_prey);
    
    int nRows = prey_rm.rows();
    int nCols = prey_rm.cols();
    
    matrix<Type> kernel_real, kernel_imag;
    buildComplexKernel(ft_pred_kernel_real, ft_pred_kernel_imag, kernel_real, kernel_imag);
    
    vector<std::complex<Type> > fft_in = prepareFFTInput(prey_rm);
    
    auto fft_forward = atomic::batched_fft_b(fft_in);
    
    auto fft_forward_scaled = scaleFFTOutput(fft_forward, nCols);
    
    auto mult = multiplyByKernel(fft_forward_scaled, kernel_real, kernel_imag, nRows, nCols);
    
    vector<Type> fft_inverse = atomic::batched_fft_b(mult, true).real();
    
    auto prey_matrix = convertFFTOutputToMatrix<Type>(fft_inverse, nRows, nCols);
    
    // Extract desired block: columns (w_full.size()-w.size()) to end.
    int start_col = static_cast<int>(w_full.size()) - static_cast<int>(w.size());
    typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixRM;
    MatrixRM ret_block = prey_matrix.block(0, start_col, nRows, w.size());
    
    ret_block = ret_block.array() * convertToRowMajor(search_vol).array();
    
    matrix<Type> fish_mask = n ;
    matrix<Type> f = feeding_level ;
    //matrix<Type> ret = ret_block;
    
    for(int j = 0 ; j < fish_mask.rows() ; j++) {
      for (int k = 0 ; k < fish_mask.cols() ; k++ ){
        fish_mask(j , k) = (n(j , k) > 0) ? 1 : 0 ;
        f(j , k) = (1 - feeding_level(j , k))*fish_mask(j , k) ;
        //ret(j , k) = f(j , k)*ret_block(j,k) ;
      }
    }
    
    ret_block = ret_block.array() * convertToRowMajor(f).array();
    
    matrix<Type> ret = ret_block;
    //return ret;
    
    for(int j = 0 ; j < ret_block.rows() ; j++) {
      for (int k = 0 ; k < ret_block.cols() ; k++ ){
        diet(j , k , i) = ret_block(j , k) ;
        //ret(j , k) = f(j , k)*ret_block(j,k) ;
      }
    }
    
  }
  return diet ;
  
}

template <class Type>
matrix <Type> propDiet(
    array<Type> diet ,
    vector<Type> w
) {
  matrix<Type> tmp(diet.rows() , diet.cols()) ; 
  tmp.setZero() ;
  matrix<Type> ret = tmp ;
  
  for (int i=0; i<ret.cols(); ++i) {
    matrix<Type> d = diet.col(i).matrix() ;
    vector<Type> tot = d.rowwise().sum() ;
    tmp.col(i) = tot ;
  }
  
  for  (int i=0; i<ret.rows(); ++i) {
    vector<Type> tot = tmp.rowwise().sum() ;
    for  (int j=0; j<ret.cols(); ++j) {
      ret(i , j) = tmp(i , j)/tot(i) ;
    }
  }
  return ret ;
}

template <class Type>
matrix <Type> propDiet2(
    array<Type> diet ,
    vector<Type> w ,
    vector<int> min_size
) {
  matrix<Type> tmp(diet.rows() , diet.cols()) ; 
  tmp.setZero() ;
  matrix<Type> ret = tmp ;
  
  for (int i=0; i<ret.cols(); ++i) {
    matrix<Type> d = diet.col(i).matrix() ;
    for (int j=0; j < ret.rows() ; j++) {
      vector<Type> v = d.row(j) ;
      //int idx = static_cast<int>(w.size()) - (min_size(j));
      tmp(j , i) = v.tail(int(w.size()) - min_size(j)).sum() ;
    }
  }
  
  for  (int i=0; i<ret.rows(); ++i) {
    vector<Type> tot = tmp.rowwise().sum() ;
    for  (int j=0; j<ret.cols(); ++j) {
      ret(i , j) = tmp(i , j)/tot(i) ;
    }
  }
  return ret ;
}

/// survey stuff - Not needed by me currently
//////////////// new bit

template<class Type>
array<Type> get_biomass_at_size(array<Type> saved_n,vector<Type> w_dw){
  array<Type> ret = saved_n;
  for (int i=0; i<ret.cols(); ++i){
    ret.col(i) = saved_n.col(i) * w_dw(i);
  }
  return(ret);
}

template<class Type>
array<Type> get_biomass_at_size_at_beg_year(array<Type> saved_n,vector<Type> w_dw, Type dt, int years){
  array<Type> ret(saved_n.rows(),saved_n.dim(1),years);
  int cur_year =1;
  Type cur_time =0;
  int time_save = 0;
  while(cur_year < (Type(years) + 1e-5)){
    cur_time += dt;
    if (cur_time >= (Type(cur_year) - 1e-5)){
      int tmp_ind=cur_year - int(1);
      ret.col(tmp_ind) = get_biomass_at_size(saved_n.col(time_save),w_dw);
      cur_year += 1;
    }
    time_save += int(1);
  }
  return(ret);
}

template<class Type>
array<Type> get_biomass_dist(array<Type> biomass_at_size , int knife_edge){
  array<Type> ret(biomass_at_size.rows(),biomass_at_size.dim(1),biomass_at_size.cols()) ;
  ret.setZero() ;
  
  array<Type> Q(biomass_at_size.rows() , biomass_at_size.dim(1)) ;
  Q.setZero() ;
  for(int i=0 ; i < Q.rows() ; ++i){
    for (int j=0 ; j < Q.cols() ; ++j){
      Q(i,j) =  (j >= knife_edge) ? 1 : 0 ;
    }
  }
  
  for(int i=0 ; i < biomass_at_size.cols() ; ++i){
    ret.col(i) = biomass_at_size.col(i) * Q ;
    vector<Type> total = ret.col(i).matrix().rowwise().sum() ;
    for(int j=0 ; j < total.size() ; ++j) {
      for (int k = 0 ; k < biomass_at_size.dim(1) ; ++k) {
        ret(j , k , i) = ret(j , k , i)/total(j) ;
      }
    }
  }
  return ret ;
}
// 
// template<class Type>
// matrix<Type> getQ(vector<Type> sizeClasses, vector<Type> knifeedge,vector<Type> qmax){
//   // Create object for output
//   int nsize = sizeClasses.size();
//   int nspecies = knifeedge.size();
//   matrix<Type> Qmatrix(nspecies,nsize);
//   Qmatrix.setZero();
//   
//   // loop over species
//   for(int i=0; i < nspecies; i++) {
//     // loop over sizes
//     for(int j=0; j < nsize; j++) {
//       
//       // condition
//       if(sizeClasses[j] >= knifeedge[i]){
//         Qmatrix(i,j) = exp(qmax[i]);
//       }
//     }
//   }
//   
//   return(Qmatrix);
// }
// 
// template<class Type>
// vector<Type> multiplyQ(matrix<Type> Qmatrix, array<Type> biomass){
//   
//   // Create object for output
//   int nspecies = Qmatrix.rows();
//   vector<Type> surveyBiomass(nspecies);
//   surveyBiomass.setZero();
//   
//   
//   // // loop over species
//   for(int i=0; i < nspecies; i++) {
//     
//     for (int j=0; j < Qmatrix.cols(); j++){
//       surveyBiomass(i) += Qmatrix(i,j) * biomass(i,j);
//     }
//   }
//   
//   return(surveyBiomass);
// }



template <class Type>
Type objective_function<Type>::operator()() {
  // Load data
  DATA_VECTOR(species_params_interaction_resource); // May consider estimating this
  // DATA_MATRIX(species_params_interaction); 
  DATA_VECTOR(w_full);
  DATA_VECTOR(w);
  DATA_VECTOR(dw_full);
  DATA_MATRIX(ft_pred_kernel_real);
  DATA_MATRIX(ft_pred_kernel_imag);
  DATA_MATRIX(search_vol);
  DATA_MATRIX(intake_max);
  DATA_VECTOR(species_params_alpha);
  //DATA_MATRIX(metab);
  DATA_VECTOR(k) ; // used to calculate metabolism (in case you want to toggle this)
  DATA_VECTOR(ks) ; // used to calculate metabolism (in case you want to toggle this)
  DATA_VECTOR(p) ; // used to calculate metabolism (in case you want to toggle this)
  DATA_MATRIX(psi);
  DATA_VECTOR(dw);
  DATA_MATRIX(ft_pred_kernel_p_real);
  DATA_MATRIX(ft_pred_kernel_p_imag);
  DATA_MATRIX(ft_mask);
  DATA_ARRAY(dummy) ; // This is catchability matrix - only works when all species fished by unique gear
  
  DATA_ARRAY(selectivity);
  //DATA_ARRAY(catchability); // Im estimating this internally
  
  DATA_MATRIX(mu_b);
  DATA_VECTOR(species_params_erepro);
  DATA_IVECTOR(w_min_idx);
  DATA_VECTOR(species_params_rmax);
  DATA_VECTOR(rr_pp);
  DATA_VECTOR(cc_pp);
  DATA_SCALAR(dt);
  
  //DATA_VECTOR(effort);  /// need to make this a matrix
  DATA_VECTOR(n_pp);
  DATA_MATRIX(n);
  
  //PARAMETER_VECTOR(log_qmax);  // number of species
  //DATA_VECTOR(knifeedge_w); 
  //DATA_ARRAY(log_survey);
  
  //For time varying external mortality
  //DATA_ARRAY(Zselectivity); // array of mu_b for each species -- estimating this now
  DATA_ARRAY(Zcatchability); // matrix of n_spec X n_spec, diag = 1
  
  DATA_INTEGER(times);
  
  /// Data Needed for Fit
  DATA_MATRIX(log_catches);
  DATA_MATRIX(log_trawl); // trawl catch
  DATA_ARRAY(Tselectivity); // trawl selectivity
  DATA_MATRIX(trawl_effort) ; //trawl effort
  
  // Data needed for sex split reproduction
  DATA_VECTOR(sex_split);
  DATA_VECTOR(female_ratio);
  DATA_IVECTOR(repro_idx);
  DATA_IVECTOR(min_size) ;
  DATA_INTEGER(knife_edge) ;
  DATA_IVECTOR(w_max_idx) ;
  
  //Fitting Info
  DATA_ARRAY(obs_dist) ;
  DATA_ARRAY(obs_diet) ;
  
  // Estimated Parameters
  // PARAMETER_VECTOR(log_resource_interaction) ; 
  // vector<Type> species_params_interaction_resource = exp(log_resource_interaction) ;
  PARAMETER_MATRIX(log_effort);
  matrix<Type> effort = exp_mat_ele(log_effort); 
  PARAMETER_MATRIX(log_Zscale);
  matrix<Type> Zscale = exp_mat_ele(log_Zscale); 
  PARAMETER_MATRIX(log_r);
  matrix<Type> Rscale = exp_mat_ele(log_r); 
  PARAMETER_ARRAY(log_Tcatchability) ;
  array<Type> Tcatchability = exp_arr_ele(log_Tcatchability);
  PARAMETER_VECTOR(log_catchability) ;
  array<Type> catchability = vec_to_diag(dummy , exp(log_catchability)) ;
  
  // Parameter SDs
  // PARAMETER_VECTOR(log_sdobs);
  // vector<Type> sdobs = exp(log_sdobs);
  PARAMETER_VECTOR(log_sdtrawl) ;
  vector<Type> sdtrawl = exp(log_sdtrawl) ;
  
  //Interaction Matrix
  PARAMETER_MATRIX(log_interaction) ;
  matrix<Type> species_params_interaction = exp_mat_ele(log_interaction);
  
  //External Mortality 
  PARAMETER_VECTOR(log_z0pre) ;
  vector<Type> z0pre = exp(log_z0pre) ;
  PARAMETER_VECTOR(log_z0exp) ;
  vector<Type> z0exp = exp(log_z0exp) ;
  array<Type> Zselectivity = set_allo_mort(w , z0pre , z0exp) ;
  
  // store n
  array<Type> n_save(n.rows(),n.cols(),times); //
  matrix<Type> new_n = n;
  vector<Type> new_npp = n_pp;
  vector<Type> ret_rm(n.rows()*n.cols() + n_pp.size());
  array<Type> catches(n.rows(),times);
  array<Type> ann_catches(n.rows(),log_effort.rows());
  array<Type> biomass(n.rows(),times);
  array<Type> ann_biomass(n.rows(),effort.rows());
  array<Type> trawl(n.rows(),times);
  array<Type> ann_trawl(n.rows(),effort.rows()) ;
  array<Type> diet_prop(n.rows(),n.rows()+1,times);
  array<Type> ann_diet(n.rows() , n.rows()+1,effort.rows()) ;
  
  // fishing mortality
  array<Type> f_mort = sort_annual_effort(selectivity, catchability, effort,times,dt);
  array<Type> z_mort = sort_annual_z(Zselectivity, Zcatchability, Zscale,times,dt);
  
  
  //for use in the ffts hopefully the names rnt important elsewhere - Need to figure out what these are
  DATA_INTEGER(nRows);
  DATA_INTEGER(nCols);
  set_fft_dims((int)nRows,(int)nCols);
  
  // Trawl Catch 
  array<Type> trawl_catch = sort_annual_effort(Tselectivity , Tcatchability , trawl_effort , times , dt) ;
  
  // Recruitment Anomaly
  matrix<Type> r_scale = sort_annual_r(Rscale,times,dt); 
  
  // useful thing
  vector <Type> w_dw = w * dw;
  
  // Toggle Size Fitting
  DATA_INTEGER(size_dist_fit) ;
  DATA_INTEGER(diet_dist_fit) ;
  
  
  for(int j=0; j< times; ++j){
    ret_rm = running_model(species_params_interaction_resource, new_npp, new_n,  
                           species_params_interaction, w_full, w, dw_full,ft_pred_kernel_real, 
                           ft_pred_kernel_imag, search_vol,intake_max,
                           species_params_alpha, //metab, 
                           psi, dw, ft_pred_kernel_p_real, ft_pred_kernel_p_imag,
                           ft_mask, f_mort.col(j), 
                           mu_b,species_params_erepro, w_min_idx, species_params_rmax, rr_pp,cc_pp,dt ,
                           k, ks,  p ,
                           z_mort.col(j) , trawl_catch.col(j) , 
                           repro_idx, sex_split, female_ratio ,
                           vector<Type>(r_scale.col(j))
                           );
    new_n = get_N_after(ret_rm,n.rows(),n.cols());
    new_npp = get_Npp_after(ret_rm,n.rows(),n.cols());
    n_save.col(j) = new_n;
    catches.col(j) = get_catch(new_n, f_mort.col(j),w_dw,dt);
    biomass.col(j) = get_biomass(new_n ,w_dw,dt) ;
    trawl.col(j) = get_catch(new_n , trawl_catch.col(j) , w_dw , dt) ;
    matrix <Type> new_encounter = getEncounter(species_params_interaction_resource, new_npp, new_n, species_params_interaction, w_full, w, dw_full,ft_pred_kernel_real, ft_pred_kernel_imag,search_vol);
    matrix <Type> new_feeding_level = FeedingLevel(intake_max, new_encounter);
    array<Type> new_diet = getDiet(species_params_interaction_resource,new_npp, new_n, species_params_interaction, w_full, w, dw_full,ft_pred_kernel_real, ft_pred_kernel_imag,search_vol , new_feeding_level);
    diet_prop.col(j) = propDiet2(new_diet , w , min_size);
  }
  
  ann_catches = get_ann_catches(catches,dt,log_effort.rows());
  ann_biomass = get_ann_biomass(biomass,dt,effort.rows());
  ann_trawl = get_ann_catches(trawl,dt,effort.rows()) ;
  ann_diet = get_ann_diet(diet_prop , dt , effort.rows()) ;
  
  ///survey
  array<Type> biomass_at_size(n.rows(),n.cols(),log_effort.rows());
  biomass_at_size = get_biomass_at_size_at_beg_year(n_save,w_dw,dt,log_effort.rows());
  array<Type> biomass_dist(n.rows(),n.cols(),log_effort.rows());
  biomass_dist = get_biomass_dist(biomass_at_size , knife_edge) ;
  
  REPORT(biomass_at_size);
  REPORT(biomass_dist) ;
  REPORT(n_save);
  REPORT(catches);
  REPORT(ann_catches);
  REPORT(f_mort);
  REPORT(z_mort) ;
  REPORT(ann_biomass) ;
  REPORT(ann_trawl) ;
  REPORT(diet_prop) ;
  REPORT(ann_diet) ;
  //REPORT(allo_mort) ;
  
  //Type sdobs= 0.1;
  Type sdZ = 0.2;
  Type sdR = 1 ;
  Type sdobs= 0.1;
  // Type sdtrawl = 0.1 ;
  
  
  // Add your model code here
  Type nll = 0.0;
  
  //int skipValue=5;
  // I need to try and find how to remove female snowcrab from this (index 6-1). I will get an error if catch = 0
  for (int i=0; i<log_catches.rows(); ++i){
    if(i == 5) continue;
    for (int j=0; j<log_catches.cols(); ++j){
      nll +=   -dnorm(log_catches(i,j),log(ann_catches(i,j)),sdobs,true); // for parameter sdobs, sdobs(i)
    }
  }
  
  for (int i=0; i<log_trawl.rows(); ++i){
    for (int j=0; j<log_trawl.cols(); ++j){
      nll +=   -dnorm(log_trawl(i,j),log(ann_trawl(i,j)),sdtrawl(i),true); // for parameter sdtrawl, sdtrawl(i)
    }
  }
  
  // should this be here?
  for (int j=0; j<log_Zscale.cols(); ++j){
    nll +=   -dnorm(log_Zscale(0,j), Type(0.0), sdZ,true);
  }
  
  for (int i=1; i<log_Zscale.rows(); ++i){
    for (int j=0; j<log_Zscale.cols(); ++j){
      nll +=   -dnorm(log_Zscale(i,j), log_Zscale(i-1,j), sdZ,true);
    }
  }
  
  for (int i = 0; i<log_r.rows();++i){
    for(int j = 0; j<log_r.cols();++j){
      nll +=   -dnorm(log_r(i,j), Type(0.0),sdR,true); 
    }
  } 
  
  
  if(diet_dist_fit == 1) {
    // diet proportion (does this need to be log transformed)
    for (int i=1; i< obs_diet.cols(); ++i) { //skipping first year because I don't have data for it
      for (int j=1; j < 3; ++j) { // hard code which species you're comparing observations for, 
        for (int k=0; k < obs_diet.dim(1); ++k){ 
          if( !R_IsNA(asDouble(obs_diet(j,k,i)))) {
            nll +=   -dnorm(obs_diet(j,k,i),ann_diet(j,k,i),sdobs,true); // for parameter sdobs, sdobs(i)
          }
        }
      }
    }
  }

  if(size_dist_fit == 1) {
    // size distribution of trawl - all species but capelin & shrimp (does this need to be log transformed)
    for (int i=0; i< obs_dist.cols(); ++i) {
      for (int j=1; j < 5; ++j) { // hard code which species you're comparing observations for, it may be wise to ultimately change this to just compare at observed sizes
        for (int k= knife_edge; k < w_max_idx(j) ; ++k){ //it may be wise to ultimately change this to just compare at observed sizes, k = knife_edge , add vector for max_w
          nll +=   -dnorm(obs_dist(j,k,i),biomass_dist(j,k,i),sdobs,true); // for parameter sdobs, sdobs(i)
        }
      }
    }
  }
  
  
  
  
  // for (int i = 0; i<log_r.rows(); ++i){
  //   for (int j=0; j < log_r.cols(); ++j){
  //     nll +=   -dnorm(log_r(i,j), Type(0.0), sdR(j),true);
  //   }
  // }
  
  // for (int j=0; j<log_r.cols(); ++j){
  //   nll +=   -dnorm(log_r(0,j), Type(0.0), sdR,true);
  // }
  // 
  // for (int i=1; i<log_r.rows(); ++i){
  //   for (int j=0; j<log_r.cols(); ++j){
  //     nll +=   -dnorm(log_r(i,j), log_r(i-1,j), sdR,true);
  //   }
  // }
  
  return nll;
  
  
  
  
  
  ///survey
  // array<Type> biomass_at_size(n.rows(),n.cols(),log_effort.rows());
  // biomass_at_size = get_biomass_at_size_at_beg_year(n_save,w_dw,dt,log_effort.rows());
  // 
  // 
  // 
  // REPORT(biomass_at_size);
  // 
  // matrix<Type> Qmatrix = getQ(w, knifeedge_w, log_qmax);
  // 
  // matrix<Type> predSurvey(knifeedge_w.size(), log_effort.rows());
  // predSurvey.setZero();
  // 
  // for (int y = 0; y < log_effort.rows(); y++) {
  //   
  //   // Extract biomass for year y 
  //   array<Type> biomass_y = biomass_at_size.col(y);
  //   
  //   //Rcout << biomass_y.dim << "\n";
  //   
  //   // Get survey predictions for year y
  //   vector<Type> predSurvey_y = multiplyQ(Qmatrix, biomass_y);
  //   predSurvey.col(y) = predSurvey_y;
  // }
  // 
  // REPORT(predSurvey);
  // REPORT(Qmatrix);
  
  //  for (int i=0; i<log_survey.rows(); ++i){
  //	  for (int j=0; j<log_survey.cols(); ++j){
  //		nll +=   -dnorm(log_survey(i,j),log(predSurvey(i,j)),sdobs,true);  /// pred survey is species by years
  //	  } 
  //  }
  
  
  //return nll;
}

