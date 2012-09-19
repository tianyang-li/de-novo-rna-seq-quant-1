// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef ACCPM_BLAS_INTERFACE_H
#define ACCPM_BLAS_INTERFACE_H

#include "AccpmVector.h"
#include "AccpmGenMatrix.h"
using namespace Accpm;

/**
 * @file AccpmBlasInterface.h
 * @ingroup AccpmLA
 */

/**
 * Functions for interfacing with C/C++ Blas/Lapack functionality.
 * Currently we only use lapack++ for all blas/lapack level functions.
 */
/**
 * Level 1 BLAS operations
 */

/** Dot Product of vectors dx and dy */
double AccpmLADotProd(const AccpmVector &dx, const AccpmVector &dy);

/** Scale vector dx by da */
void AccpmLAScale(double da, AccpmVector &dx);

/** Vector scaling: dy = da * dx */
void AccpmLAMult(AccpmVector &dy, double da, const AccpmVector &dx);

/**  
 * Combined vector scaling and addition:  dy = dy + da * dx 
 * saxpy
*/
void AccpmLAAddMult(AccpmVector &dy, double da, const AccpmVector &dx);

/** \brief 1-Norm
 *
 * Returns the sum of the absolute values: \f$|x|_1=\sum_i|x_i|\f$
 */
double AccpmLANorm1(const AccpmVector &dx);

/** \brief 2-Norm, Euclidean Norm
 *
 * Returns the euclidean norm of the vector:
 * \f$|x|_2=\sqrt{\sum_i|x_i|^2}\f$
 */
double AccpmLANorm2(const AccpmVector &dx);

/** \brief Infinity-Norm
 *
 * Returns the Infinity norm of a vector, which is the absolute value
 * of its maximum element: \f$|x|_{\infty}=\max_i|x_i|\f$
 */
double AccpmLANormInf(const AccpmVector &x);

/**
 * Level 2 BLAS operations
 */

/** Perform the matrix-vector operation y := alpha*A'*x + beta*y */
void AccpmLAMatTransVecMult(const RealMatrix &A, 
			    const AccpmVector &dx, 
			    AccpmVector &dy,
			    double alpha = 1.0, double beta = 0.0);

/** Perform the matrix-vector operation y := alpha*A*x + beta*y */
void AccpmLAMatVecMult(const RealMatrix &A, 
		       const AccpmVector &dx, 
		       AccpmVector &dy, 
		       double alpha = 1.0, double beta = 0.0);

/** Perform the rank 1 operation A := alpha*dx*dy' + A */
void AccpmLAR1Update(RealMatrix &A, const AccpmVector &dx, 
		     const AccpmVector &dy, double alpha = 1.0);

void AccpmLAR1Update(SymmetricMatrix &A, const AccpmVector &dx,
		    double alpha = 1.0);
/**
 * Level 3 BLAS operations
 */
/** Perform the matrix-matrix operation C := alpha*A*B + beta*C */
void AccpmLAMatMatMult(const RealMatrix &A, 
		       const RealMatrix &B, RealMatrix &C, 
		       double alpha = 1.0, double beta = 0.0)
;
/** Perform the matrix-matrix operation C := alpha*A'*B + beta*C */
void AccpmLAMatTransMatMult(const RealMatrix &A, 
			    const RealMatrix &B, RealMatrix &C, 

			    double alpha = 1.0, double beta = 0.0);
/** Perform the matrix-matrix operation C := alpha*A*B' + beta*C */
void AccpmLAMatMatTransMult(const RealMatrix &A, 
			    const RealMatrix &B, RealMatrix &C, 
			    double alpha = 1.0, double beta = 0.0);

/** Perform  the  symmetric  rank  k  operations
 * C := alpha*A*A’ + beta*C
 * Only the lower traingular part of C is used
 */
void AccpmLAR1Update(SymmetricMatrix &C, RealMatrix &A,
		      double alpha = 1.0, double beta = 1.0);
#endif
