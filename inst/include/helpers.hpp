// Helpers: https://kaskr.github.io/adcomp/matrix_arrays_8cpp-example.html

// =========================== convertDimsSeals ===============================
template<class Type>
vector <Type> convertDimsUpNpp(vector<Type> n_pp, vector<Type> w , vector<Type> w_new){
  vector<Type> ret(w_new.size());
  ret.setZero();
  vector<Type> tmp = n_pp.tail(w.size());
  for (int i=0; i<tmp.size(); ++i){
    ret(i) = tmp(i);
  }
  return ret;
}

extern "C" SEXP call_convertDimsUpNpp(SEXP n_pp , SEXP w , SEXP w_new){
  vector <double > y = convertDimsUpNpp (asVector<double>(n_pp) , asVector<double>(w) , asVector<double>(w_new));
  return asSEXP (y);
}

// =========================== convertDimsUpNpp ===============================
template<class Type>
matrix <Type> convertDimsUpN(matrix<Type> n , vector<Type> w_new){
  matrix<Type> ret(n.rows() , w_new.size());
  ret.setZero();
  ret.block(0, 0 , ret.rows(), n.cols()) += n;
  return ret;
}

extern "C" SEXP call_convertDimsUpN(SEXP n , SEXP w_new){
  matrix <double > y = convertDimsUpN (asMatrix<double>(n) , asVector<double>(w_new));
  return asSEXP (y);
}

// =========================== convertDimsDownNpp ===============================
template<class Type>
vector <Type> convertDimsDownNpp(vector<Type> n_pp, vector<Type> w , vector<Type> npp_convert){
  vector<Type> ret(n_pp.size());
  ret.setZero();
  ret.tail(w.size()) += npp_convert.head(w.size());
  return ret;
}

extern "C" SEXP call_convertDimsDownNpp(SEXP n_pp , SEXP w , SEXP npp_convert){
  vector <double > y = convertDimsDownNpp (asVector<double>(n_pp) , asVector<double>(w) , asVector<double>(npp_convert));
  return asSEXP (y);
}

// =========================== convertDimsDownN ===============================
template<class Type>
matrix <Type> convertDimsDownN(matrix<Type> n , matrix<Type> n_convert){
  matrix<Type> ret(n.rows() , n.cols());
  ret.setZero();
  ret = n_convert.block(0, 0 , ret.rows(), ret.cols());
  return ret;
}

extern "C" SEXP call_convertDimsDownN(SEXP n , SEXP n_convert){
  matrix <double > y = convertDimsDownN (asMatrix<double>(n) , asMatrix<double>(n_convert));
  return asSEXP (y);
}

//// Other useful functions
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

// make everything in an array an exponential for 3D array
template<class Type>
array<Type> exp_arr3_ele(array<Type> M_in){
  array<Type> ret = M_in;
  for (int i = 0; i < M_in.rows(); i++) {
    for (int j = 0; j < M_in.cols(); j++) {
      for (int k = 0; k < M_in.dim(1); k++) {
        ret(i,k,j) = exp(M_in(i,k,j)) ;
      }
    }
  }
  return(ret);
}

