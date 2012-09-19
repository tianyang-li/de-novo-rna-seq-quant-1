// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "MCFSmoothOracleFunction.h"
#include "AccpmBlasInterface.h"

MCFSmoothOracleFunction::MCFSmoothOracleFunction(const MCFData *data) : _data(data)
{
}

MCFSmoothOracleFunction::~MCFSmoothOracleFunction()
{
}

int 
MCFSmoothOracleFunction::eval(const AccpmVector &y, AccpmVector &functionValue, 
			      AccpmGenMatrix &subGradients, AccpmGenMatrix *info)
{
  if (_data->_type == MCFData::KLEINROCK) {
    //f2 = - sum(2*sqrt(y .* C1) - 1 - y .* C1); % objective function
    int n = y.size();
    AccpmVector yc3 = _data->_capacity;
    yc3.times(y);
    AccpmVector tmp(n);
    tmp = 1;
    AccpmLAAddMult(tmp, 1, yc3);
    tmp.negate();
    for (int i = 0; i < n; ++i) {
      tmp(i) += 2*sqrt(yc3(i));
    }
    functionValue = -tmp.sum();
    
    //df2   =   - (sqrt(C1./ y ) - C1); % computes the  gradiant 
    AccpmVector c3 = _data->_capacity;
    tmp = c3;
    tmp.rdivide(y);
    for (int i = 0; i < n; ++i) {
      tmp(i) = sqrt(tmp(i));
    }
    AccpmLAAddMult(c3, -1, tmp);
    memcpy(subGradients.addr(), c3.addr(), sizeof(double)*n);

    if (info) {
      assert(info->size(0) == n);
      //d2f2  =   0.5 * sqrt(C1) .* y .^(-1.5) ;% computes the diag of the hessian
      tmp.rdivide(y);
      AccpmLAScale(0.5, tmp);
      memcpy(info->addr(), tmp.addr(), sizeof(double)*n);
    }
  } 
  
  return 0;
}
