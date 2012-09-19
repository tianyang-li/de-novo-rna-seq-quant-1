// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "MCFOracle.h"
#include "AccpmDefs.h"
#include "AccpmBlasInterface.h"

MCFOracle::MCFOracle(OracleFunction *f1, OracleFunction *f2, const MCFData *data)
  : Oracle(f1, f2), _data(data)
{
  _lowerBound = ACCPM_MINUS_INF;
}

MCFOracle::~MCFOracle()
{
}

bool
MCFOracle::computesBound() const
{
  return true;
}

double
MCFOracle::computeLowerBound() 
{
  double lowerBound = ACCPM_MINUS_INF;
  const AccpmVector &capacity = _data->_capacity;
  const AccpmVector &sol = _data->_sol;
  if (sol.size() > 0) {
    bool updateBound = true;
    for (int i = 0; i < capacity.size(); ++i) {
      //if (DBL_LT(capacity(i), sol(i), 1e-8)) {
      if (DBL_LT(capacity(i), sol(i))) {
	updateBound = false;
	break;
      }
    }
    if (updateBound) {
      if (_data->_type == MCFData::LINEAR) {
	const AccpmVector &cost = _data->_cost;
	lowerBound = -AccpmLADotProd(cost, sol);
	//_data->_stop = 1;
      } else if (_data->_type == MCFData::KLEINROCK) {
	//LowerBound = - sum(DataS.D.sol ./ (C1(:,3) -  DataS.D.sol));	
	AccpmVector tmp1 = _data->_capacity;
	AccpmLAAddMult(tmp1, -1, sol);
	AccpmVector tmp2 = sol;
	tmp2.rdivide(tmp1);
	lowerBound = -tmp2.sum();
      }
    }
  }
  if (_lowerBound < lowerBound) {
    _lowerBound = lowerBound;
  } else {
    lowerBound =  _lowerBound;
  }
  return lowerBound;
}

double
MCFOracle::getObjectiveFunctionBound() const
{
  return _lowerBound;
}
