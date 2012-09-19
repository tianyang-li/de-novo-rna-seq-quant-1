// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "AccpmLASolve.h"
#include "lapackd.h"
#include "laslv.h"

namespace Accpm {
 
void
AccpmLALinearSolve(const RealMatrix &A, RealMatrix &x, const RealMatrix &b)
{
  LaLinearSolve(A, x, b);
}

int
AccpmLASymmLinSolve(const AccpmGenMatrix &A, RealMatrix &X, const RealMatrix &B)
{
  char uplo = 'L';
  long int n = A.size(0);
  long int lda = A.inc(0) * A.gdim(0);
  long int nrhs = B.size(1);
  long int ldb = B.inc(0) * B.gdim(0);
  long int info;
  AccpmGenMatrix L(A);
  X.inject(B);

  F77NAME(dposv) (&uplo, &n, &nrhs, &L(0,0), &lda, &X(0,0), &ldb, &info);

  if (info > 0) {
    std::cerr << "AccpmLASymmLinSolve: Matrix is not positive-definite." 
	      << std::endl;
  } else if (info < 0) {
    std::cerr << "AccpmLASymmLinSolve: argument " << -info << " is invalid" 
	      << std::endl;
  }

  return info;
}

int
AccpmLASymmLinSolve(SymmetricMatrix &A, RealMatrix &X, const RealMatrix &B)
{
  LaLowerTriangMatDouble L(A);
  char uplo = 'L';
  long int n = A.size(0);
  long int lda = A.inc(0) * A.gdim(0);
  long int nrhs = B.size(1);
  long int ldb = B.inc(0) * B.gdim(0);
  long int info;
  
  X.inject(B);

  F77NAME(dposv) (&uplo, &n, &nrhs, &L(0,0), &lda, &X(0,0), &ldb, &info);
 
  if (info > 0) {
    std::cerr << "AccpmLASymmLinSolve: Symmetric Matrix is not positive-definite." 
	      << std::endl;
  } else if (info < 0) {
    std::cerr << "AccpmLASymmLinSolve: argument " << -info << " is invalid" 
	      << std::endl;
  }
  return info;
}

int
AccpmLACholeskyFactor(const RealMatrix &A, RealMatrix &L)
{
  L.copy(A);

  char uplo = 'L';
  long int n = A.size(0);
  long int lda = A.inc(0) * A.gdim(0);
  long int info;
  F77NAME(dpotrf) (&uplo, &n, &L(0,0), &lda, &info);
  if (info > 0) {
    std::cerr << "AccpmLACholeskyFactor: Matrix is not Symmetric-positive-definite." 
	      << std::endl;
  } else if (info < 0) {
    std::cerr << "AccpmLACholeskyFactor: argument " << -info << " is invalid" 
	      << std::endl;
  }
  return info;
}

void
AccpmLALinSolve(const RealMatrix &A, bool cholesky, RealMatrix &X, const RealMatrix &B)
{
  if (!cholesky) {
    AccpmLALinearSolve(A, X, B);
    return;
  }
  char uplo = 'L';
  long int n = A.size(0);
  long int lda = A.inc(0) * A.gdim(0);
  long int nrhs = B.size(1);
  long int ldb = B.inc(0) * B.gdim(0);
  long int info;
  
  X.inject(B);
  //double alpha = 1;
  //F77NAME(dtrsm)("L","U","T","N",&n,&nrhs,&alpha,&A(0,0),&n,&X(0,0),&n);
  //F77NAME(dtrsm)("L","U","N","N",&n,&nrhs,&alpha,&A(0,0),&n,&X(0,0),&n);
 
  F77NAME(dpotrs)(&uplo, &n, &nrhs, &A(0,0), &lda, &X(0,0), &ldb, &info);
  if (info < 0) {
    std::cerr << "AccpmLALinSolve: argument " << -info << " is invalid" 
	      << std::endl;
  }
  
}

}
