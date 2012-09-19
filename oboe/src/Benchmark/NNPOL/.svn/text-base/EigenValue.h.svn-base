// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef EIGENVALUE_H
#define EIGENVALUE_H
/**
 * @file EigenValue.h
 * Support for computing Eigen Values of a complex Hermitian Matrix.
 */

#define LA_COMPLEX_SUPPORT
#include "lacomplex.h"
#include "gmc.h"
#include "lavc.h"
#include "AccpmVector.h"
typedef LaGenMatComplex AccpmComplexMatrix;
typedef LaVectorComplex AccpmComplexVector;

extern "C"
{
  void F77NAME(zheev)(char *jobz, char *uplo, integer *n, doublecomplex *a, integer *lda, double *w, 
		       doublecomplex *work, integer *lwork, double *rwork, integer *info);
}

int computeEigenValues(const AccpmComplexMatrix &A, Accpm::AccpmVector &eigenValues);

#endif
