// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "AccpmGenMatrix.h"
#include <iostream>
#include "AccpmBlasInterface.h"
#include "blaspp.h"

using std::cout;
using std::endl;

namespace Accpm 
{
  
AccpmGenMatrix::AccpmGenMatrix() : RealMatrix()
{
}
    
AccpmGenMatrix::AccpmGenMatrix(int m, int n)
  : RealMatrix(m,n)
{
}
    
AccpmGenMatrix::AccpmGenMatrix(double*v, int m, int n, bool row_ordering) 
  : RealMatrix(v, m, n, row_ordering)
{
  assert(row_ordering == false); // We only support column ordered matrices for now
}

AccpmGenMatrix::AccpmGenMatrix(const AccpmGenMatrix &rhs) 
  : RealMatrix(rhs)
{
}

AccpmGenMatrix::AccpmGenMatrix(const RealMatrix &rhs) 
  : RealMatrix(rhs)
{
}

AccpmGenMatrix::AccpmGenMatrix(const AccpmVector &rhs) 
  : RealMatrix(rhs)
{
}
/*
AccpmGenMatrix 
AccpmGenMatrix::operator()(const LaIndex& I, const LaIndex& J) const
{
  return RealMatrix::operator()(I,J);
}
    
AccpmGenMatrix 
AccpmGenMatrix::operator()(const LaIndex& I, const LaIndex& J) 
{
  return RealMatrix::operator()(I,J);
}
*/
AccpmGenMatrix*
AccpmGenMatrix::transpose() const
{
  AccpmGenMatrix *transposeMat = new AccpmGenMatrix(size(1), size(0));
  for (int i = 0; i < size(1); ++i) {
    transposeMat->assignRow(i, getColumn(i));
  }
  return transposeMat;
}

void 
AccpmGenMatrix::assignColumn(int colId, const AccpmVector &v)
{
  int numCols = size(1);
  assert(colId >= 0);
  assert(colId < numCols);
  assert(size(0) == v.size());
  double *destaddr = addr() + size(0)*(colId); 
  memcpy(destaddr, v.addr(), sizeof(double)*v.size());
}

void 
AccpmGenMatrix::assignRow(int rowId, const AccpmVector &v)
{
  assert(size(1) == v.size());
  for (int j = 0; j < size(1); ++j) {
    (*this)(rowId, j) = v(j);
  }
}

AccpmGenMatrix& 
AccpmGenMatrix::operator=(double s)
{

  RealMatrix::operator=(s);
  return *this;
}

AccpmGenMatrix& 
AccpmGenMatrix::operator=(const AccpmGenMatrix &s)
{
  RealMatrix::operator=(s);
  return *this;
}

RealVector
AccpmGenMatrix::getColumn(int i) const
{ 
  
  LaIndex colIndex(i,i);
  LaIndex rowIndex = index(0);
  AccpmVector col = RealMatrix::operator()(rowIndex, colIndex);
  return col;
}

RealVector
AccpmGenMatrix::getRow(int i) const
{
  LaIndex rowIndex(i,i);
  LaIndex colIndex = index(1);
  AccpmVector row = RealMatrix::operator()(rowIndex, colIndex);
  return row;
}

void
AccpmGenMatrix::scaleColumn(int i, double d) 
{
  /*LaIndex colIndex(i,i);
  LaIndex rowIndex = index(0);
  AccpmVector col = (RealMatrix::operator()(rowIndex, colIndex));
  AccpmLAScale(d, col);
  */
  integer n = size(0);
  integer incx = 1;
  
  F77NAME(dscal)(&n, &d, addr() + (i * n), &incx);
}

void
AccpmGenMatrix::scaleColumn(int i, const AccpmVector &d) 
{
  assert(size(0) == d.size());
#if (defined(WIN32))
  for (int ix = 0; ix < size(0); ++ix) {
	(*this)(ix, i) *= d(ix);
  }
#else
  LaIndex colIndex(i,i);
  LaIndex rowIndex = index(0);
  AccpmVector col(RealMatrix::operator()(rowIndex, colIndex));
  col.times(d);
#endif
}

void 
AccpmGenMatrix::addMult(double scale, const AccpmGenMatrix &b)
{
  assert(size(0) == b.size(0));
  assert(size(1) == b.size(1));
  integer n = size(0)*size(1);
  integer inc = 1;
  F77NAME(daxpy)(&n, &scale, b.addr(), &inc, addr(), &inc);
}

void
AccpmGenMatrix::scale(double scale)
{
  integer n = size(0)*size(1);
  integer inc = 1;
  
  F77NAME(dscal)(&n, &scale, addr(), &inc);
  
}

}
