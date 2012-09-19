// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "AccpmVector.h"
#include "AccpmDefs.h"

#include "AccpmBlasInterface.h"

#include <iostream>

namespace Accpm {

AccpmVector::AccpmVector() : RealVector()
{
}

AccpmVector::AccpmVector(int n) : RealVector(n)
{
}

AccpmVector::AccpmVector(double *v, int n) : RealVector(v, n)
{
}

AccpmVector::AccpmVector(const LaGenMatDouble &s) : RealVector(s)
{
}

AccpmVector&
AccpmVector::operator=(double s)
{
  RealVector::operator=(s);
  return *this;
}

AccpmVector& 
AccpmVector::operator=(const LaGenMatDouble &s)
{
  RealVector::operator=(s);
  return *this;
}

bool 
AccpmVector::operator==(const AccpmVector &v) const
{
  if (size() != v.size()) {
    return false;
  }
  for (int i = 0; i < size(); ++i) {
    if (!DBL_CMP((*this)(i), v(i))) {
      return false;
    }
  }
  return true;
}

void 
AccpmVector::copy(const StdRealVector &v)
{
  for (unsigned int i = 0; i < v.size(); ++i) {
    (*this)(i) = v[i];
  }
}

void 
AccpmVector::append(const AccpmVector &v)
{
  AccpmVector tmp(size() + v.size());
  if (tmp.size() > 0) {
    memcpy(tmp.addr(), addr(), sizeof(double) * size());
    memcpy(tmp.addr() + size(), v.addr(), sizeof(double) * v.size());
  }
  ref(tmp);
}

void 
AccpmVector::append(const double element)
{
  int length = size();
  AccpmVector tmp(length + 1);
  if (length > 0) {
    memcpy(tmp.addr(), addr(), sizeof(double) * length);
  }
  tmp(length) = element;
  ref(tmp);
}

int
AccpmVector::deleteElem(int id)
{
  int length = size();
  assert(0 <= id && id < length);
  AccpmVector tmp(length - 1);
  if (length > 0) {
    memcpy(tmp.addr(), addr(), sizeof(double) * id);
    memcpy(tmp.addr()+id, addr()+id+1, sizeof(double) * (length-1-id));
  }
  ref(tmp);
  return id;
}

void
AccpmVector::negate()
{
  AccpmLAScale(-1, *this);
}

void 
AccpmVector::invert()
{
  for (int i = 0; i < size(); ++i) {
    double val = (*this)(i);
    assert(!DBL_CMP(val, 0));
    (*this)(i) = 1.0/val;
  }
}

void 
AccpmVector::times(const AccpmVector &v)
{
  assert(size() == v.size());
  for (int i = 0; i < size(); ++i) {
    (*this)(i) *= v(i);
  }
}

void 
AccpmVector::rdivide(const AccpmVector &v)
{
  assert(size() == v.size());
  for (int i = 0; i < size(); ++i) {
    double val = v(i);
    assert(!DBL_CMP(val, 0));
    (*this)(i) /= val;
  }
}

double
AccpmVector::sum() const
{
  AccpmVector ones(size());
  ones = 1;
  return AccpmLADotProd(*this, ones);
}

double
AccpmVector::min(int *index) const
{
  double returnValue = ACCPM_PLUS_INF;
  for (int i = 0; i < size(); ++i) {
    if (returnValue > (*this)(i)) {
      returnValue = (*this)(i);
      if (index) {
	*index = i;
      }
    }
  }
  return returnValue;
}

double
AccpmVector::max(int *index) const
{
  double returnValue = ACCPM_MINUS_INF;
  for (int i = 0; i < size(); ++i) {
    if (returnValue < (*this)(i)) {
      returnValue = (*this)(i);
      if (index) {
	*index = i;
      }
    }
  }
  return returnValue;
}

}


