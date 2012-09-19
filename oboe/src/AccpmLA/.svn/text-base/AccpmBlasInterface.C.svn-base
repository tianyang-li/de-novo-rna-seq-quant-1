// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team

#include "AccpmBlasInterface.h"
#include "blaspp.h"

/** Calls Lapack++ Blas_Dot_Prod */
double 
AccpmLADotProd(const AccpmVector &dx, const AccpmVector &dy)
{
  return Blas_Dot_Prod(dx, dy);
}

/** Calls Lapack++ Blas_Scale */
void 
AccpmLAScale(double da, AccpmVector &dx)
{
  Blas_Scale(da, dx);
}

void
AccpmLAMult(AccpmVector &dy, double da, const AccpmVector &dx)
{
  Blas_Mult(dy, da, dx);
}

void 
AccpmLAAddMult(AccpmVector &dy, double da, const AccpmVector &dx)
{
  Blas_Add_Mult(dy, da, dx);
}

double 
AccpmLANorm1(const AccpmVector &dx)
{
  return Blas_Norm1(dx);
}

double 
AccpmLANorm2(const AccpmVector &dx)
{
  return Blas_Norm2(dx);
}

void 
AccpmLAMatTransVecMult(const RealMatrix &A, 
		       const AccpmVector &dx, 
		       AccpmVector &dy,
		       double alpha, double beta)
{
  Blas_Mat_Trans_Vec_Mult(A, dx, dy, alpha, beta);
}

void 
AccpmLAMatVecMult(const RealMatrix &A, 
		  const AccpmVector &dx, 
		  AccpmVector &dy, 
		  double alpha, double beta)
{
  Blas_Mat_Vec_Mult(A, dx, dy, alpha, beta);
}

void 
AccpmLAR1Update(RealMatrix &A, const AccpmVector &dx, 
		const AccpmVector &dy, double alpha)
{
  Blas_R1_Update(A, dx, dy, alpha);
}

#include <cmath>
double 
AccpmLANormInf(const AccpmVector &x)
{   
  integer index = Blas_Index_Max(x);
  return std::fabs(x(index));
}

void 
AccpmLAR1Update(SymmetricMatrix &A, const AccpmVector &dx,
		double alpha)
{
  Blas_R1_Update(A, dx, alpha);
}

void 
AccpmLAMatMatMult(const RealMatrix &A, 
		  const RealMatrix &B, RealMatrix &C, 
		  double alpha, double beta)
{
  Blas_Mat_Mat_Mult(A, B, C, alpha, beta);
}

void
AccpmLAMatTransMatMult(const RealMatrix &A, 
		       const RealMatrix &B, RealMatrix &C, 
		       double alpha, double beta)
{
  Blas_Mat_Trans_Mat_Mult(A, B, C, alpha, beta);
}

void 
AccpmLAMatMatTransMult(const RealMatrix &A, 
		       const RealMatrix &B, RealMatrix &C, 
		       double alpha, double beta)
{
  Blas_Mat_Mat_Trans_Mult(A, B, C, alpha, beta);
}

void
AccpmLAR1Update(SymmetricMatrix &C, RealMatrix &A,
		double alpha, double beta)
{
  Blas_R1_Update(C, A, alpha, beta);
}
