// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef MCF_ORACLE_H
#define MCF_ORACLE_H

#include "MCFOracleFunction.h"

class MCFOracle : public Oracle {

 private:
  const MCFData *_data;
  double _lowerBound;

 public:
  MCFOracle(OracleFunction *f1, OracleFunction *f2, const MCFData *data);
  virtual ~MCFOracle();
  virtual bool computesBound() const;
  virtual double getObjectiveFunctionBound() const;
  virtual double computeLowerBound();
};

#endif
