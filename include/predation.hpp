//// functions involving FFT
template <typename T>
inline double get_base(const T &x) {
  return x;
}
template <typename T>
inline double get_base(const CppAD::AD<T>& x) {
  return get_base( CppAD::Value(x) );
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

template<class T>
vector<std::complex<T> > cplx(vector<T> x) {
  vector<std::complex<T> > xc(x.size());
  for (int i = 0; i < x.size(); i++) {
    xc[i] = std::complex<T>(x[i], T(0));
  }
  return xc;
}

// =========================== computePrey ===============================
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

extern "C" SEXP call_computePrey(SEXP species_params_interaction_resource, SEXP n_pp , SEXP n , SEXP species_params_interaction , SEXP w_full , SEXP w , SEXP dw_full){
  matrix <double > y = computePrey (asVector<double>(species_params_interaction_resource),asVector <double>(n_pp) , asMatrix<double>(n),asMatrix<double>(species_params_interaction),asVector<double>(w_full),asVector<double>(dw_full) );
  return asSEXP (y);
}

// =========================== computeQMatrix ===============================
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

extern "C" SEXP call_computeQMatrix(SEXP n, SEXP w_full , SEXP w , SEXP feeding_level , SEXP search_vol , SEXP dw){
  matrix <double > y = computeQMatrix (asMatrix<double>(n),asVector <double>(w_full),asVector<double>(w),asMatrix<double>(feeding_level),asMatrix<double>(search_vol),asVector<double>(dw) );
  return asSEXP (y);
}

// =========================== getPredRate ===============================
template <class Type>
matrix<Type> getPredRate(
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

extern "C" SEXP call_getPredRate(SEXP n, SEXP w_full , SEXP w , SEXP feeding_level , SEXP search_vol , SEXP dw , SEXP ft_pred_kernel_p_real , SEXP ft_pred_kernel_p_imag , SEXP ft_mask){
  matrix <double > y = getPredRate (asMatrix<double>(n),asVector <double>(w_full),asVector<double>(w),asMatrix<double>(feeding_level),asMatrix<double>(search_vol),asVector<double>(dw),asMatrix<double>(ft_pred_kernel_p_real),asMatrix<double>(ft_pred_kernel_p_imag),asMatrix<double>(ft_mask) );
  return asSEXP (y);
}

// =========================== getEncounter ===============================
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

  // Multiply elementwise by search_vol.
  ret_block = ret_block.array() * convertToRowMajor(search_vol).array();

  matrix<Type> ret = ret_block;
  return ret;
}

extern "C" SEXP call_getEncounter(SEXP n, SEXP w_full , SEXP w , SEXP feeding_level , SEXP search_vol , SEXP dw){
  matrix <double > y = getEncounter (asMatrix<double>(n),asVector <double>(w_full),asVector<double>(w),asMatrix<double>(feeding_level),asMatrix<double>(search_vol),asVector<double>(dw),asMatrix<double>(ft_pred_kernel_p_real),asMatrix<double>(ft_pred_kernel_p_imag),asMatrix<double>(search_vol) );
  return asSEXP (y);
}

// =========================== getFeedingLevel ===============================
template<class Type>
matrix <Type> getFeedingLevel(matrix<Type> intake_max, matrix<Type> encounter){
  matrix<Type> ret(encounter.rows(),encounter.cols()) ;
  for (int i=0; i<ret.rows(); ++i){
    for (int j=0; j<ret.cols(); ++j){
      ret(i,j) = encounter(i,j) / (encounter(i,j) + intake_max(i,j));
    }
  }
  return ret;
}

// =========================== getPredMort ===============================
template<class Type>
matrix <Type> PredMort(vector<Type> w_full, vector<Type> w, matrix<Type> interaction, matrix<Type> pred_rate){
  matrix<Type> ret(pred_rate.rows(),w.size()) ;
  matrix <Type> tmp= pred_rate.block(0,w_full.size()-w.size(),pred_rate.rows(),w.size());
  ret = interaction.transpose() * tmp; // ask Gustav about drop==FALSE
  return ret;
}

// =========================== getResourceMort ===============================
template<class Type>
vector <Type> ResourceMort(vector<Type> species_params_interaction_resource, matrix<Type> pred_rate){
  vector<Type> ret(pred_rate.cols()) ;
  matrix<Type> tmp = pred_rate.transpose();
  ret = tmp * species_params_interaction_resource;
  return ret;
}
