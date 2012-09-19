// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "QPNDHSmoothOF.h"
#include "AccpmBlasInterface.h"

int 
QPNDHSmoothOF::eval(const AccpmVector &y, AccpmVector &functionValue, 
		    AccpmGenMatrix &subGradients, AccpmGenMatrix *info) 
{

  assert(_nx + _ny == y.size());
  AccpmVector x = y(LaIndex(0,_nx-1));
  
  functionValue = 0.5*AccpmLADotProd(x,x);
  subGradients = 0;
  memcpy(subGradients.addr(), y.addr(), sizeof(double)*_nx);
  if (info) {
    assert(info->size(0) == y.size());
    *info = 0;
    for(int i = 0; i < _nx; ++i) {
      (*info)(i,0) = 1;
    }
  }
  
  return 0;
}

 
