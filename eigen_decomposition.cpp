// [[Rcpp::depends(lidR)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <SpatialIndex.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace lidR;

//' Fast Eigenvalues decomposition for k nearest neighbors using a C++ function
//'
//' C++ helper function to compute eigenvalues for geometric feature
//' calculation.
//'
//' @param las LAS element
//' @param k k nearest neighbors
//' @param ncpu number of cpu cores to use
//' @return The function returns for every point the 3 eigenvalues and the
//' third element of the third eigenvector. These values are needed to compute
//' planarity, linerity, verticality etc. in the add_geometry function
//' @author Julian Frey <julian.frey@@iww.uni-freiburg.de>
//' @seealso \link{add_geometry}
//' @export eigen_decomposition
// [[Rcpp::export]]
NumericMatrix eigen_decomposition(S4 las, int k, int ncpu = 1)
{
  DataFrame data = as<DataFrame>(las.slot("data"));
  NumericVector X = data["X"];
  NumericVector Y = data["Y"];
  NumericVector Z = data["Z"];
  unsigned int npoints = X.size();

  NumericMatrix out(npoints, 4);

  SpatialIndex index(las);

#pragma omp parallel for num_threads(ncpu)
  for (unsigned int i = 0 ; i < npoints ; i++)
  {
    arma::mat A(k,3);
    arma::mat coeff;  // Principle component matrix
    arma::mat score;
    arma::vec latent; // Eigenvalues in descending order

    PointXYZ p(X[i], Y[i], Z[i]);

    std::vector<PointXYZ> pts;
    index.knn(p, k, pts);

    for (unsigned int j = 0 ; j < pts.size() ; j++)
    {
      A(j,0) = pts[j].x;
      A(j,1) = pts[j].y;
      A(j,2) = pts[j].z;
    }

    arma::princomp(coeff, score, latent, A);

#pragma omp critical
{
  out(i, 0) = latent[0];
  out(i, 1) = latent[1];
  out(i, 2) = latent[2];
  out(i, 3) = coeff[8];
}
  }

  return out;
}
