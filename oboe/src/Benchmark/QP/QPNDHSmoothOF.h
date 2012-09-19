// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef QPNDHSMOOTH_OF_H
#define QPNDHSMOOTH_OF_H

#include "Oracle.h"

using namespace Accpm;

class QPNDHSmoothOF : public OracleFunction {

private:
  int _nx;
  int _ny;

public:
  QPNDHSmoothOF(int nx, int ny) : OracleFunction(), _nx(nx), _ny(ny){};
  

  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info);

};

#endif
