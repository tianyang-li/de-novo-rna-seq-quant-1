// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef ACCPM_DEFS_H
#define ACCPM_DEFS_H

#include <float.h>
#define ACCPM_PLUS_INF DBL_MAX/100
#define ACCPM_MINUS_INF -ACCPM_PLUS_INF
#include <iostream>
#include <cmath>

#include <set>
typedef std::set<int> StdIntSet;

inline bool IS_FINITE(double x) { return (x < ACCPM_PLUS_INF && x > ACCPM_MINUS_INF); }
//inline bool DBL_CMP(double x, double y) { return x == y; }
//inline bool DBL_LT(double x, double y) { return (y - x > 0); }

inline bool 
DBL_CMP(double A, double B, double precision = 1e-15)
{
  if (A == B)
    return true;

  if (A == 0) {
    return std::fabs(B) < 1e-30;
  }
  if (B == 0) {
    return std::fabs(A) < 1e-30;
  }

  double relativeError;
  if (std::fabs(B) > std::fabs(A))
    relativeError = std::fabs((A - B) / B);
  else
    relativeError = std::fabs((A - B) / A);
  if (relativeError <= precision)
    return true;
  return false;
}

inline bool 
DBL_LT(double A, double B, double precision = 1e-14)
{
  if (A == B || DBL_CMP(A, B, precision))
    return false;
  return A < B;
}

inline bool 
DBL_GT(double A, double B, double precision = 1e-14)
{
  return (DBL_LT(B, A, precision));
}


struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

inline void AccpmWarning(const char *msg) { std::cout << "WARNING:" << msg << std::endl; }
inline void AccpmError(const char *msg) { std::cout << "ERROR:" << msg << std::endl; }

#endif
