// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

//[[Rcpp::export]]
Eigen::VectorXd expmAtv_cpp(Eigen::MatrixXd evecs, Eigen::VectorXd evals,
  Eigen::VectorXd v, Eigen::VectorXd d, Eigen::VectorXd dInv,
  double t = 300, bool preMultiply = false) {

    MatrixXd unscaled = evecs *
      (evals.array() * t).exp().matrix().asDiagonal() *
      evecs.transpose();

    // correct small values
    double* raw = unscaled.data();
    double n = unscaled.rows() * unscaled.cols();
    for(int i=0; i<n; i++) {
      double tmp = *raw;
      *(raw++) = std::abs(tmp);
    }

    if(preMultiply) {
      return v.transpose() * d.asDiagonal() * unscaled * dInv.asDiagonal();
    } else {
      return d.asDiagonal() * unscaled * dInv.asDiagonal() * v;
    }
}

//[[Rcpp::export]]
Rcpp::List expm_cpp(Eigen::MatrixXd A, double delta = 1e-15,
                    double t = 300, double tol = 1e-15) {

  int N = A.rows();

  VectorXd diag = A.diagonal();
  VectorXd dsuper = A.diagonal(1);
  VectorXd dsub = A.diagonal(-1);

  // make A symmetrizable by adding small noise
  if(delta > 0) {

    // shift main diagonal
    double delta2 = 2 * delta;
    for(int i=0; i < N; i++)
      diag(i) -= delta2;

    // shift super and sub diagonals
    for(int i=0; i<N-1; i++) {
      dsuper(i) += delta;
      dsub(i) += delta;
    }
    dsuper(0) += delta;
    dsub(N-2) += delta;

  }

  // direct construction of similarity-transformed symmetric version of matrix
  VectorXd E(N-1);
  for(int i=0; i<N-1; i++) {
    E(i) = sqrt(dsuper(i) * dsub(i));
  }

  // LAPACK internals
  char compz = 'I';
  MatrixXd evecs(N,N);
  VectorXd work(2*N-2);
  int info = 0;

  F77_CALL(dsteqr)(&compz, &N, diag.data(), E.data(), evecs.data(), &N,
                   work.data(), &info);

  // // LAPACK internals
  // char jobz = 'V';
  // char range = 'A';
  // double dzero  = 0.0;
  // int izero = 0;
  // double abstol = tol;
  // VectorXi isuppz(2 * N);
  // double work_query;
  // int info = 0, lwork = -1, liwork = -1, iwork_query;
  // int M = N;
  //
  // // output
  // VectorXd evals(N);
  // MatrixXd evecs(N,N);
  //
  // // workspace query
  // F77_CALL(dstevr)(&jobz, &range, &N, diag.data(), E.data(), &dzero, &dzero,
  //                  &izero, &izero, &abstol, &N, evals.data(), evecs.data(), &N,
  //                  isuppz.data(), &work_query, &lwork, &iwork_query, &liwork,
  //                  &info);
  //
  // lwork = (int) work_query;
  // liwork = (int) iwork_query;
  // VectorXd work(lwork);
  // VectorXi iwork(liwork);
  //
  // // eigen decomposition
  // F77_CALL(dstevr)(&jobz, &range, &N, diag.data(), E.data(), &dzero, &dzero,
  //                  &izero, &izero, &abstol, &M, evals.data(), evecs.data(), &N,
  //                  isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
  //                  &info);
  //

   // build similarity transform
   VectorXd d(N), dInv(N);
   d(0) = 1.0;
   dInv(0) = 1.0;
   for(int i=1; i<N; i++) {
     double tmp = std::sqrt(dsub(i-1) / dsuper(i-1)) * d(i-1);
     d(i) = tmp;
     dInv(i) = 1/tmp;
   }

   MatrixXd unscaled = evecs *
     (diag.array() * t).exp().matrix().asDiagonal() *
     evecs.transpose();

   // "correct" small, negative values
   double* raw = unscaled.data();
   int entries = N*N;
   for(int i=0; i<entries; i++) {
     double tmp = *raw;
     *(raw++) = std::abs(tmp);
   }

  MatrixXd expm = d.asDiagonal() * unscaled * dInv.asDiagonal();

  return Rcpp::List::create(
    Rcpp::Named("expm") = expm,
    Rcpp::Named("vectors") = evecs,
    Rcpp::Named("values") = diag,
    Rcpp::Named("d") = d,
    Rcpp::Named("dInv") = dInv
  );
}
