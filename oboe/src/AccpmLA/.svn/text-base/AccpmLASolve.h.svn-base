// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef ACCPMLA_SOLVE_H
#define ACCPMLA_SOLVE_H

#include "AccpmVector.h"
#include "AccpmGenMatrix.h"

namespace Accpm {

/**
 * Solves Ax = b using Lapack++ interface to Lapack.
 **/
void AccpmLALinearSolve(const RealMatrix &A, RealMatrix &x, const RealMatrix &b);

/**
 * Solves AX = B.
 **/
int AccpmLASymmLinSolve(const AccpmGenMatrix &A, RealMatrix &X, const RealMatrix &B);

/**
 * Solves AX = B, for Symmetric Matrices.
 **/
int AccpmLASymmLinSolve(SymmetricMatrix &A, RealMatrix &X, const RealMatrix &B);

/**
 * Computes the Cholesky factor for matrix A and returns it in L.
 * @returns The Cholesky factorization L of A 0 LL^T
 * @returns 0 on success, or
 * info value from dpotrf.
 */
int AccpmLACholeskyFactor(const RealMatrix &A, RealMatrix &L);

/**
* Solves AX = B. It uses the Cholesky factorization of A computed by a
* previous call to AccmpLACholeskyFactor.
*/
void AccpmLALinSolve(const RealMatrix &A, bool cholesky, RealMatrix &X, const RealMatrix &B);

}
#endif
