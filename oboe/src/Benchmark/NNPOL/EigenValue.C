// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "EigenValue.h"
#include "symd.h"
#include "laslv.h"


using namespace Accpm;

void
test()
{
  int n = 5;
  AccpmComplexMatrix A(n,n);
  A = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      LaComplex c(j-i+1,0);
      A(i,j) = c;
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      A(i,j) = A(j,i);
    }
  }
  std::cout << "A:\n" << A << std::endl;
  AccpmVector eigenValues(n);

  if (!computeEigenValues(A, eigenValues)) {
    std::cout << "The Eigen Values are:\n" << eigenValues << std::endl;
  }
  LaSymmMatDouble S(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      S(i,j) = j-i+1;
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      S(i,j) = S(j,i);
    }
  }
  std::cout << "S:\n" << S << std::endl;
  AccpmVector eig(n);
  LaEigSolve(S, eig);
  std::cout << "The Eigen Values are:\n" << eig << std::endl;
}

int
computeEigenValues(const AccpmComplexMatrix &A, AccpmVector &eigenValues)
{
  char uplo = 'U';
  integer n = A.size(0);
  assert(A.size(0) == A.size(1));
  AccpmComplexMatrix copyA(A);
  AccpmVector D(n);
  AccpmVector E(n-1);
  AccpmVector tau(n-1);
  integer info;
  integer lwork = (40+2)*n;
  AccpmComplexMatrix work(lwork,1);
  AccpmVector rwork(3*n-2);
  char compz = 'N';
 
  F77NAME(zheev)(&compz, &uplo, &n, &copyA(0,0), &n, &D(0), &work(0,0), &lwork, &rwork(0), &info);
  if (info == 0) {
    //std::cout << "zheev successful." << std::endl;
    //std::cout << "D:\n" << D << std::endl;
    eigenValues = D;
  } else if (info < 0) {
    std::cout << "zheev unsuccessful argument " << -info << " has illegal value." << std::endl;
  } else {
    std::cout << "zheev unsuccessful : The algorithm did not converge." << std::endl;
  }
  
  return info;
}

