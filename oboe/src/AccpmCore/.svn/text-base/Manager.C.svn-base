// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//
//

#include "Manager.h"
#include "Oracle.h"
#include "AccpmDynMatrix.h"
#include "AccpmBlasInterface.h"

#include <cmath>

using std::endl;

namespace Accpm {

/**
 * PointGen Constructor. Called when Manager is created.
 * Creates Q vector with appropriate dimensions, 
 * but ATQA is created as an empty matrix, it gets initialized when needed.
 */
PointGen::PointGen(const Parameters &param) : _q(param.getIntParameter("NumVariables")), _ATQA()
{
  _ws = _zs = _ss = 1;
  _rho = param.getRealParameter("Rho");
} 

/**
 * Basic intitialization routine used by the Constructor
 */
void
Manager::init()
{
  _phase2 = false;
  _prevPhase2 = false;
  _currentPointIsFeasible = false;
  _numOuterIteration = 0;
  _numInnerIteration = 0;
  _exitCode = ITERATING;
  _rho = _param->getRealParameter("Rho");
  _objectiveFunction =  ACCPM_PLUS_INF;
  _objLB = ACCPM_MINUS_INF;
  _objUB = ACCPM_PLUS_INF;
  _relativeGap = ACCPM_PLUS_INF;
  _weightEpigraphCut = _param->getRealParameter("WeightEpigraphCutInit");
  _f2 = 0;
  _df2 = 0;
  _d2f2 = 0;
  if (getB()) {
    _b = new AccpmVector(_param->getIntParameter("NumVariables"));
    _b->RealVector::copy(*getB());
    if (_param->getOptType() == OPT_MAX) {
      _b->negate();
    }
    _df2 = new AccpmVector(_b->size());
    _df2->RealVector::copy(*_b);
  } else {
    _b = 0;
  }
  _sizeOfs = 0;
  _currentCutId = -1;
  _numCuts = 0;
  _hasSmoothOracle = false;
}

/**
 * The Constructor for the Manager. It does some basic initialization and creates
 * vectors of appropriate dimensions for:
 * Current Point Y
 * Best Point
 * Current Point Z
 *
 * Also initializes PointGen object
 */
Manager::Manager(const Parameters *param) : _param(param),
					    _currentPointy(_param->getIntParameter("NumVariables")),
					    _bestPointy(_param->getIntParameter("NumVariables")),
					    _currentPointz(_param->getIntParameter("NumSubProblems")),
					    _pointGen(*_param)
{
  init();

 
  //_activeCuts(false, true, 0, 0);
  
  if (_param->getStartingPoint()) {
    _currentPointy = *_param->getStartingPoint();
  } else {
    _currentPointy = 0;
  }

  _bestPointy = 0;
  _currentPointz = 0;
  _objLB = _param->getRealParameter("ObjectiveLB");
  _objUB = _param->getRealParameter("ObjectiveUB");
 
  if (IS_FINITE(_objUB)) {
    updateRelativeGap();
  }

  if (_param->getOracle()) {
    _hasSmoothOracle = _param->getOracle()->hasSmoothOracle();
  }
  updateQ();
  // _pointGen._q = _rho* _param->getIntParameter("Proximal");
}

Manager::~Manager() {
 
  delete _b;
  delete _df2;
  delete _d2f2;
    
  for (CutSet::iterator iter = _cutSet.begin(); iter != _cutSet.end();) {
    const AccpmVectorIntPair *v = (*iter).first;
    ++iter;
    delete v->first; // delete the AccpmVector
	delete v;
  }
}

const AccpmVector& 
Manager::getCurrentY() const
{ 
  return _currentPointy; 
}

const AccpmVector& 
Manager::getActiveY() const
{
  return _currentPointy; 
}

const AccpmVector& 
Manager::getBestY(bool activeVar) const
{
  return _bestPointy;
}

int 
Manager::getActiveCuts(AccpmGenMatrix &cuts, bool activeVar) const 
{ 
  cuts = _activeCuts.getM();
  return 0;
}

int 
Manager::getActiveCuts(const AccpmGenMatrix *&cuts, bool activeVar) const 
{ 
  cuts = &_activeCutsM;
  return 0;
}

/**
 * Temporary debug routine for some sanity checking fof Manager. To be removed.
 */
void
Manager::check()
{
  int n = _bestPointy.size();
  if (n > 0) {
    _phase2 = true;
    int numCuts = 6;
    AccpmVector v(n);
    v = 0;
    for (int i = 0; i < numCuts; ++i) {
      v(1) = 1;
      _activeCuts.addColumn(v);
    }
    AccpmVector v1(numCuts);
    v1(0) = 1;
    v1(1) = 2;
    v1(2) = 3;
    v1(3) = 0;
    v1(4) = 1;
    v1(5) = 0;
    _subProblemIndex = v1;
    v1 = 1;
    _rhsCoef = v1;
    _cutOccurrence.append(v1);
    std::cout << "ActiveCuts:\n" << _activeCuts.getM()  << std::endl;
    _pointGen._q = 1*_rho* _param->getIntParameter("Proximal");
  }
}

/**
 * Update the Manager with the new cut information for the non-smooth Oracle f1
 **/
int 
Manager::update1(const AccpmVector &functionValue, const AccpmGenMatrix &subGradients, 
		 const AccpmGenMatrix &subProblemIndex, const AccpmVector *upperFunctionValues)
{
  ++_numOuterIteration;
  _exitCode = ITERATING;
  int numCuts = subGradients.size(1);
  if (numCuts) {
    assert(numCuts == functionValue.size());
    assert(numCuts == subProblemIndex.size(0));
  }
  assert(subProblemIndex.size(1) == 1);
  AccpmVector index(subProblemIndex.size(0));
  index.RealVector::inject(subProblemIndex);
  _currentPointIsFeasible = index.sum() > 0;
  
  if (_currentPointIsFeasible) {
    _phase2 = _currentPointIsFeasible;
    updateQ();
  }

  OptType optType = _param->getOptType();

  AccpmVector functionVal;
  functionVal = functionValue;
  AccpmVector upperFunctionVal;

  if (upperFunctionValues) {
    upperFunctionVal = *upperFunctionValues;
  } else {
    upperFunctionVal = functionValue;
  }

  AccpmGenMatrix newCuts(subGradients);
  if (_currentPointIsFeasible && optType == OPT_MAX) {
    functionVal.negate();
    upperFunctionVal.negate();
    
    for (int i = 0; i < newCuts.size(1); ++i) {
      newCuts.scaleColumn(i, -1);
    }
  }
  int numCutsFiltered = processCuts(functionVal, newCuts, index, upperFunctionVal);
   
  if (!_hasSmoothOracle) {
    update();
  }
    
  return numCutsFiltered;
}

/**
 * Add the new cut information from the non-smooth Oracle f1.
 * It is called by update1() function which takes care of transforming
 * the cut information for a minimization problem.
 * Update the internal structures of Manager corresponding to the new information.
 * Returns number of filtered cuts.
 **/
int 
Manager::processCuts(const AccpmVector &functionValue, const AccpmGenMatrix &newCuts, 
		     const AccpmVector &subProblemIndex, const AccpmVector &upperFunctionValues)
{
  int numNewCuts = newCuts.size(1);
  
  double correction = 0;
  if (numNewCuts && _phase2) {
    if (_param->getIntParameter("ConvexityCheck") || _param->getIntParameter("ConvexityFix")) {
      correction = convexityFix(functionValue, newCuts, subProblemIndex, upperFunctionValues);
    }
  }
 
  // Update objective function
  if (_currentPointIsFeasible) {
    _objectiveFunction = AccpmLADotProd(functionValue, (*_param->getPi()));
    if (_b) {
      _objectiveFunction += AccpmLADotProd(_currentPointy, *_b);
    }
    _objectiveFunction += correction;
  }
  //update(_currentPointIsFeasible);
  int numCutsFiltered = 0;
  // Add the cuts with or without filtering
  if (numNewCuts) {
    for (int i = 0; i < numNewCuts; ++i) {
	  AccpmVector *cut = new AccpmVector(newCuts.getColumn(i));
      AccpmVectorIntPair *cutp = new AccpmVectorIntPair(cut, (int)subProblemIndex(i));
 
      if (_param->getIntParameter("Filter")) {
	    int cutId = findCut(*cutp);
	    if (cutId < 0) {
	     double newRhsCoef = - functionValue(i) + AccpmLADotProd(*cut, _currentPointy);
	     addCut(cutp, (int)subProblemIndex(i), 1, newRhsCoef);
	    } else {
	      numCutsFiltered++;
	     //std::cout << "Filtered Cut: " << i << "\n" << *cut << std::endl;
	     //std::cout << "Original Cut: " << cutId << "\n" << _activeCuts.getM().getColumn(cutId) << std::endl;
		  double newRhsCoef = - functionValue(i) + AccpmLADotProd(*cut, _currentPointy);
	      if (DBL_GT(_rhsCoef(cutId), newRhsCoef)) {
	       _rhsCoef(cutId) = newRhsCoef;
	      }
	      _cutOccurrence(cutId)++;
		  delete cut;
	      delete cutp;	
	    }
      } else {
	_activeCuts.addColumn(*cut);
	++_numCuts;
	_subProblemIndex.append(subProblemIndex(i));
	_cutOccurrence.append(1);
	double norm = AccpmLANorm2(*cut);
	if (subProblemIndex(i) > 0) {
	  ++norm;
	}
	_activeCutNorm.append(norm);
	delete cut;
	delete cutp;
      }
    }
    
    _activeCutsM = _activeCuts.getM();
    
    if (_param->getIntParameter("Verbosity") > 1 && numCutsFiltered) {
      std::cout << numCutsFiltered << " of the " << numNewCuts << " new cuts were filtered" << std::endl;
    }
    
    // Update rhs if we added new cuts
    if (!_param->getIntParameter("Filter")) {
      updateRhs(functionValue, newCuts);
    }
  }

  return numCutsFiltered;
}

int
Manager::addCut(const AccpmVectorIntPair *cutp, int subProblemIndex, double weight, double rhs)
{ 
  ++_currentCutId;
  if (_param->getIntParameter("Filter")) {
    _cutSet[cutp] = _currentCutId;
  }
  int idInMatrix = _activeCuts.addColumn(*(cutp->first));
  ++_numCuts;
  _subProblemIndex.append(subProblemIndex);
  _cutOccurrence.append(weight);
  double norm = AccpmLANorm2(*(cutp->first));
  if (subProblemIndex > 0) {
    ++norm;
  }
  _activeCutNorm.append(norm);
  _rhsCoef.append(rhs);
  
  return idInMatrix;
}

/**
 * Fix loose optimality cuts.
 * Applies to oracle that deliver loose optimality cuts. (E.g., if the oracle
 * computes the value f(y) by a maximization process and that process is not
 * conducted to its end). 
 *
 * Updates the Objective Upper Bound: _objUB
 *
 * Returns a value by which the objective function needs to be corrected due
 * to convexity violation.
 *
 * Must be called only in phase 2.
 */
double
Manager::convexityFix(const AccpmVector &functionValue,const AccpmGenMatrix &subGradients, 
		      const AccpmVector &subProblemIndex, const AccpmVector &upperFunctionValues) 
{

  assert(_phase2);
  double correction = 0;
  AccpmVector dy = _bestPointy;
  AccpmLAAddMult(dy, -1, _currentPointy);
  AccpmVector delta(functionValue.size());
  AccpmLAMatTransVecMult(subGradients, dy, delta, 1, 0);
  AccpmLAAddMult(delta, 1, functionValue);
  
  double objective = AccpmLADotProd(*_param->getPi(), delta);
  if (_b) {
    objective += AccpmLADotProd(_bestPointy, *_b);
  }
  if (_hasSmoothOracle) {
    AccpmWarning("The convexity analysis with smooth oracle is not yet supported");
  }
  double test2 = _objUB - objective;
  double test = test2 / std::max(1., fabs(_objUB));
  if (DBL_LT(test, -0.01 * _param->getRealParameter("Tolerance"))) {
    if (_param->getIntParameter("ConvexityCheck")) {
      _exitCode = CONVEXITY_FAILURE;
      if (_param->getOptType() == OPT_MIN) {
	AccpmError("\nThe best upper bound is underestimated by ");
      } else {
	AccpmError("\nThe best lower bound is overestimated by ");
      }
      std::cout << -test2 << std::endl;
    } else {
      _objUB = objective + 0.05 * _param->getRealParameter("Tolerance");
    }
  }

  // Check whether former cuts invalidate the point.
  for (int i = 0; i < _param->getIntParameter("NumSubProblems"); ++i) {
    AccpmVector cut;
    delta = AccpmVector();
    for (int j = 0; j < _subProblemIndex.size(); ++j) {
      if (_subProblemIndex(j) == i) {
	cut = _activeCutsM.getColumn(j);
	double del = AccpmLADotProd(cut, _currentPointy);
	del -= _rhsCoef(j);
	delta.append(del);
      }
    }
    double deltaMax = delta.max();
    test2 = -deltaMax;
    int index;
    for (index = 0; index < subProblemIndex.size(); ++index) {
      if (subProblemIndex(index) == i) {
	test2 += upperFunctionValues(index);
	break;
      }
    }
    test = test2 / std::max(1., fabs(_objUB));
    if (DBL_LT(test, -0.01 * _param->getRealParameter("Tolerance"))) {
      AccpmError("");
      if (_param->getOptType() == OPT_MIN) {
	std::cout << "The function value with index: " << i 
		  << " is underestimated by: " << -test2 << std::endl;
      } else {
	std::cout << "The function value with index: " << i 
		  << " is overestimated by: " << -test2 << std::endl;
      }
      if (_param->getIntParameter("ConvexityCheck")) {
	_exitCode = CONVEXITY_FAILURE;
      } else {
	double value = deltaMax;
	value = std::max(upperFunctionValues(index), deltaMax) + 0.1 * _param->getRealParameter("Tolerance");
	correction += upperFunctionValues(index) - functionValue(index);
      }
    }
  }
  
  return correction;
}

/**
 * Update the right hand side corresponding to the new cuts just returned by the Oracle
 */
void
Manager::updateRhs(const AccpmVector &functionValue, const AccpmGenMatrix &subGradients)
{ 
  AccpmVector rhs = functionValue;
  AccpmLAMatTransVecMult(subGradients, _currentPointy, rhs, 1, -1);
  _rhsCoef.append(rhs);
}

/**
 * Update the internal data of Manager:
 *
 * Rho
 * Phase2 flag
 * Upper Bound on Objective Function
 * Best Point
 * Relative Gap
 * Weight on the Epigraph Cut
 **/
void
Manager::update()
{
  
  if (_currentPointIsFeasible) {
    updateRho();
    if (!DBL_GT(_objectiveFunction, _objUB)) {
      _objUB = _objectiveFunction;
      _bestPointy = _currentPointy;
    }
  }
  
  if (_phase2) {
    updateRelativeGap();
  } else {
    if (!_param->getIntParameter("FixedProximalCenter")) {
      _bestPointy = _currentPointy; // only update best point
    }
  }

  updateEpigraphWeight();
}

void
Manager::updateRelativeGap()
{
  _relativeGap = (_objUB - _objLB)/std::max(1., fabs(_objUB));
  if (!DBL_GT(_relativeGap, _param->getRealParameter("Tolerance"))) {
    _exitCode = RELATIVE_GAP_REACHED;
  }
}

/**
 * Update the Proximal term.
 *
 */
void
Manager::updateRho()
{
  if (_phase2 && _param->getIntParameter("Proximal")) {
    if (_param->getIntParameter("DynamicRho")) {	
      double rhoFactor = 1;
      double expansion = 0.5;
      if (DBL_LT(_objectiveFunction, _objUB)) {
        rhoFactor = expansion;
      } else {
        rhoFactor = 1./expansion;
      }
      double newRho = _rho * rhoFactor;
      newRho = std::min(_param->getRealParameter("RhoMax"), newRho);
      newRho = std::max(_param->getRealParameter("RhoMin"), newRho);
      _rho = newRho;
      updateQ();
    } else { // allow use to change rho between iterations
      _rho = _param->getRealParameter("Rho");
      updateQ();
    }
  }
}

void
Manager::updateQ()
{
  _pointGen._q = _rho * _param->getIntParameter("Proximal");
}

void
Manager::updateEpigraphWeight()
{
  if (_phase2 && IS_FINITE(_objUB)) {
    _weightEpigraphCut += _param->getRealParameter("WeightEpigraphCutInc");
  }
}

/**
 * Update the smooth function f2 values
 */
int 
Manager::update2(double functionValue, const AccpmGenMatrix &subGradients, 
		 const AccpmGenMatrix &hessian)
{
  OptType optType = _param->getOptType();
  int n = _param->getIntParameter("NumVariables");
  assert(subGradients.size(1) == 1);

  if (!_df2) {
    _df2 = new AccpmVector(n);
  }  
  *_df2 = 0;
  
  if (optType == OPT_MAX) {
    _f2 = -functionValue;
    AccpmLAAddMult(*_df2, -1, subGradients);
  }  else {
    _f2 = functionValue;
    AccpmLAAddMult(*_df2, 1, subGradients);
  }
  
  if (!_d2f2) {
    if (_param->getIntParameter("DiagHessian")) {
      _d2f2 = new AccpmGenMatrix(n, 1);
    } else {
      _d2f2 = new AccpmGenMatrix(n, n);
    }
  }
  _d2f2->copy(hessian);
  
  if (optType == OPT_MAX) {
    for (int i = 0; i < _d2f2->size(1); ++i) {
      _d2f2->scaleColumn(i, -1);
    }
  }
  
  if (_currentPointIsFeasible) {
    _objectiveFunction += _f2;
  }

  update();

  return 0;
}


/**
 * Returns the diagonal part of the Ball Constraints
 * Mirrors InsertBallConstraint(...,'Diag')
 **/
double 
Manager::computeBallConstraint(const AccpmVector &y) const
{
  double result = 0;
  if (_param->getIntParameter("Ball")) {
    AccpmVector v = y;
    //v = y - centerBall
    AccpmLAAddMult(v, -1, *getCenterBall(true));
    double norm = AccpmLANorm2(v);
    norm *= norm;
    double radius = getRadiusBall();
    norm = std::pow(radius, 2) - norm;
    assert(norm > 0);
    result = 2.0/ norm;
  }
  return result;
}

/**
 * Returns the Rang 1  part of the Ball Constraints
 * Mirrors InsertBallConstraint(...,'Rang')
 **/
void 
Manager::computeBallConstraint(const AccpmVector &y, AccpmVector &result) const
{
  result = 0;
  if (_param->getIntParameter("Ball")) {
    result = y;
    //v = y - centerBall
    AccpmLAAddMult(result, -1, *getCenterBall(true));
    
    double result1 = computeBallConstraint(y);
    AccpmLAScale(result1, result);
  }
}

void 
Manager::computeBox1Constraint(const AccpmVector &y, AccpmVector &result) const
{
  result = 0;
  const AccpmVector *varB = getVariableLB(true);
  if (varB) {
    AccpmVector tmp = *varB;
    AccpmLAAddMult(tmp, -1, y);
    tmp.invert();
    AccpmLAAddMult(result, 1, tmp);
  }
  varB = getVariableUB(true);
  if (varB) {
    AccpmVector tmp = *varB;
    AccpmLAAddMult(tmp, -1, y);
    tmp.invert();
    AccpmLAAddMult(result, 1, tmp);
  }
}

void 
Manager::computeBox2Constraint(const AccpmVector &y, AccpmVector &result) const
{
  result = 0;
  const AccpmVector *varB = getVariableLB(true);
  if (varB) {
    AccpmVector tmp = y;
    AccpmLAAddMult(tmp, -1, *varB);
    tmp.times(tmp);
    tmp.invert();
    AccpmLAAddMult(result, 1, tmp);
  }
  varB = getVariableUB(true);
  if (varB) {
    AccpmVector tmp = *varB;
    AccpmLAAddMult(tmp, -1, y);
    tmp.times(tmp);
    tmp.invert();
    AccpmLAAddMult(result, 1, tmp);
  }
}

void
Manager::updateY(const AccpmVector &y)
{
  _currentPointy = y;
}

void 
Manager::updateVariables(const AccpmVector &y, const AccpmVector &z, const AccpmVector &s,
			 double zs, double s0, double ss,
			 const AccpmVector &sOld, double sOld0,
			 const AccpmVector &ds, double ds0)

{
  updateY(y);
  _currentPointz = z;
  _pointGen._s = s;
  _pointGen._zs = zs;
  _pointGen._ss = ss;
  _pointGen._sOld = sOld;
  _pointGen._ds = ds;
  _pointGen._x = _cutOccurrence;
  _pointGen._x.rdivide(s);
  if (_phase2) {
    _pointGen._s.append(s0);
    _pointGen._sOld.append(sOld0);
    _pointGen._ds.append(ds0);
    _pointGen._x.append(_weightEpigraphCut / s0);
  }
}

void
Manager::updateLB(double lBound)
{
  if (_param->getIntParameter("Verbosity") > 3) {
    if (DBL_LT(_objLB, lBound)) {
      std::cerr << "Lower Bound improved" << std::endl;
    } else if (!DBL_CMP(_objLB, lBound)) {
      std::cerr << "Lower Bound became worse" << std::endl;
    }
  }
  _objLB = std::max(lBound, _param->getRealParameter("ObjectiveLB"));
  updateRelativeGap();
}

int 
Manager::findCut(const AccpmVectorIntPair &v) const
{
  CutSet::const_iterator it = _cutSet.find(&v);
  if (it != _cutSet.end()) {
    return (*it).second;
  }
  return -1;
}

double
Manager::getRelativeGap() const
{
  return _relativeGap;
}

void
Manager::output(std::ostream &os) const 
{
  os << "\n-------------------------------------------------------------------------" << endl;
  os << "Exit Code:\t\t\t" << _exitCode << endl;
  os << "Optimization Type:\t\t" << _param->getOptimizationType() << endl;
  os << "Number of Variables:\t\t" << _currentPointy.size() << endl;
  os << "Number of Cuts:\t\t\t" << _numCuts << endl;
  os << "Number of Subproblems:\t\t" << _currentPointz.size() << endl;
  os << "Optimal Objective Value:\t" << _param->getOptType() * _objUB << endl;
  if (_relativeGap == ACCPM_PLUS_INF) {
    os << "Relative Gap:\t\t\tInf" << endl;
  } else {
    os << "Relative Gap:\t\t\t" << _relativeGap << endl;
  }
  os << "Number of Inner Iterations:\t" << _numInnerIteration << std::endl;
  if (_param->getIntParameter("Verbosity") > 5) {
    std::cout << "\nSolution y:\n" << _bestPointy << std::endl;
    std::cout << "\nRhs :\n" << _rhsCoef << std::endl;
  }
}
 
void 
Manager::updateActiveVariables()
{
}

void 
Manager::updateActiveVariablesForF2()
{
}

const AccpmVector*
Manager::getVariableLB(bool activeVar) const
{
  return _param->getVariableLB();
}

const AccpmVector* 
Manager::getVariableUB(bool activeVar) const
{
  return _param->getVariableUB();
}

const AccpmVector*
Manager::getB(bool activeVar) const
{
  return _param->getB();
} 

const AccpmVector*
Manager::getCenterBall(bool activeVar) const
{
  return _param->getCenterBall();
} 

double 
Manager::getRadiusBall() const
{
  return _param->getRealParameter("RadiusBall");
}

const AccpmVector *
Manager::getQ(bool activeVar) const
{
  return &_pointGen._q;
}

double 
Manager::getSmoothObj(bool activeVar) const 
{ 
  return _f2; 
}

const AccpmVector *
Manager::getSmoothGradient(bool activeVar) const 
{ 
  return _df2; 
}
    
const AccpmGenMatrix *
Manager::getSmoothHessian(bool activeVar) const 
{ 
  return _d2f2; 
}
     
int
Manager::computeSmoothComponent(bool hessian)
{
  if (_currentPointIsFeasible) {
    if (_hasSmoothOracle) {
      assert(_b == 0);
      if (hessian) {
	callSmoothOracle(_currentPointy, _f2, *_df2, _d2f2);
      } else {
	callSmoothOracle(_currentPointy, _f2, *_df2, 0);
      }
    } else {
      if (_df2) {
	_f2 = AccpmLADotProd(_currentPointy, *_df2);
	assert(_d2f2 == 0);
      }
    }
  }
  return 0;
}

int
Manager::callSmoothOracle(const AccpmVector &y, double &f2, AccpmVector &df2, AccpmGenMatrix *d2f2) const
{
  int returnValue = 0;
  OracleFunction *f = _param->getOracle()->getF2();
  if (f) {
    AccpmVector val2(1);
    AccpmGenMatrix subGrad2(df2.size(), 1);
    returnValue = f->eval(y, val2, subGrad2, d2f2);
    f2 = val2(0);
    df2.RealVector::copy(subGrad2);
    OptType optType = _param->getOptType();
    if (optType == OPT_MAX) {
      f2 *= -1;
      df2.negate();
      if (d2f2) {
	for (int i = 0; i < d2f2->size(1); ++i) {
	  d2f2->scaleColumn(i, -1);
	}
      }
    }
   
  }
  return returnValue;
}

int
Manager::getNumCuts() const
{
  return _numCuts;
}

}
