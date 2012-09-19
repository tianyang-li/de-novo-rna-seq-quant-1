// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef MCF_SMOOTH_ORACLE_FUNCTION_H
#define MCF_SMOOTH_ORACLE_FUNCTION_H

#include "MCFOracleFunction.h"

class MCFSmoothOracleFunction : public OracleFunction {
 
 private:
  const MCFData *_data;
 
 public:
  MCFSmoothOracleFunction(const MCFData *data);
  virtual ~MCFSmoothOracleFunction();
  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info);
};

#endif
