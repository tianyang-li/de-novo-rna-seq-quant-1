// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "Parameters.h"
#include "AccpmDefs.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <functional>
#include <algorithm>

using std::endl;
using std::string;
using std::cerr;
using std::cout;

namespace Accpm {

void
Parameters::init()
{
  _problemName = "OBOE General Problem"; 
  _problemAcronym = "OBOE_Problem";
  _optimizationType = "Min";
  _type = OPT_MIN;
  _numVariables = 0;
  _pi = 0;
  setNumSubProblems(1);
  
  _objLB = ACCPM_MINUS_INF;
  _objUB = ACCPM_PLUS_INF;
  
   _startingPoint = 0;
   _varLB = 0;
   _varUB = 0;
   _b = 0;
   _delta = 5;
   _eta = 0.99;
   _gamma = 0.99;
   _epsilonReal = 1e-10;
   _epsilonTol = 1e-10;

   _methodName = "ProximalACDual";
   _maxOuterIterations = 1000;
   _maxInnerIterations = 50;
   _verbosity = 1;

   _tolerance = 1e-6;
   _filter = true;

   _convexityCheck = false;
   _convexityFix = false;
   _fixedProximalCenter = false;

    _proximal = true;
    _rho = _proximal;
    _dynamicRho = false;
    _rhoMax = 100;
    _rhoMin = 1e-6;
    _weightEpigraphCutInit = 1;
    _weightEpigraphCutInc = 1;
   
    _ball = false;
    _radiusBall = 1e5;
    _centerBall = 0;
    
    _oracle = 0;
    _diagHessian = true;

    _computeLB = true;

    _checkLocSet = false;
    _lpSolverName = "GLPK";

    _avIndex = 0;
    
    _D = 0;
    _d = 0;

    setFunctionMaps();
}

void
Parameters::setFunctionMaps()
{
  // The integer valued functions
  _sintFunctionMap["NumVariables"] = &Parameters::setNumVariables;
  _sintFunctionMap["NumSubProblems"] = &Parameters::setNumSubProblems;
  _sintFunctionMap["MaxOuterIterations"] = &Parameters::setMaxOuterIterations;
  _sintFunctionMap["MaxInnerIterations"] = &Parameters::setMaxInnerIterations;
  _sintFunctionMap["Verbosity"] = &Parameters::setVerbosity;
  
  _gintFunctionMap["NumVariables"] = &Parameters::getNumVariables;  
  _gintFunctionMap["NumSubProblems"] = &Parameters::getNumSubProblems;
  _gintFunctionMap["MaxOuterIterations"] = &Parameters::getMaxOuterIterations;
  _gintFunctionMap["MaxInnerIterations"] = &Parameters::getMaxInnerIterations;
  _gintFunctionMap["Verbosity"] = &Parameters::getVerbosity;

  // The boolean functions
  _sintFunctionMap["Filter"] = &Parameters::setFilter;
  _sintFunctionMap["ConvexityCheck"] = &Parameters::setConvexityCheck;
  _sintFunctionMap["ConvexityFix"] = &Parameters::setConvexityFix;
  _sintFunctionMap["FixedProximalCenter"] = &Parameters::setFixedProximalCenter;
  _sintFunctionMap["Proximal"] = &Parameters::setProximal;
  _sintFunctionMap["DynamicRho"] = &Parameters::setDynamicRho;
  _sintFunctionMap["Ball"] = &Parameters::setBall;
  _sintFunctionMap["DiagHessian"] = &Parameters::setDiagHessian;
  _sintFunctionMap["ComputeLowerBound"] = &Parameters::setComputeLB;
  _sintFunctionMap["CheckLocSetInterior"] = &Parameters::setCheckLocSetInterior;

  _gintFunctionMap["Filter"] = &Parameters::getFilter;
  _gintFunctionMap["ConvexityCheck"] = &Parameters::getConvexityCheck;
  _gintFunctionMap["ConvexityFix"] = &Parameters::getConvexityFix;
  _gintFunctionMap["FixedProximalCenter"] = &Parameters::getFixedProximalCenter;
  _gintFunctionMap["DynamicRho"] = &Parameters::getDynamicRho;
  _gintFunctionMap["Proximal"] = &Parameters::getProximal;
  _gintFunctionMap["Ball"] = &Parameters::getBall;
  _gintFunctionMap["Box"] = &Parameters::getBox;
  _gintFunctionMap["DiagHessian"] = &Parameters::getDiagHessian;
  _gintFunctionMap["ComputeLowerBound"] = &Parameters::getComputeLB;
  _gintFunctionMap["CheckLocSetInterior"] = &Parameters::getCheckLocSetInterior;
 
  // The real valued functions
  _srealFunctionMap["ObjectiveLB"] = &Parameters::setObjLB;
  _srealFunctionMap["ObjectiveUB"] = &Parameters::setObjUB;
  _srealFunctionMap["Delta"] = &Parameters::setDelta;
  _srealFunctionMap["Eta"] = &Parameters::setEta;
  _srealFunctionMap["Gamma"] = &Parameters::setGamma;
  _srealFunctionMap["Tolerance"] = &Parameters::setTolerance;
  _srealFunctionMap["Rho"] = &Parameters::setRho;
  _srealFunctionMap["RhoMax"] = &Parameters::setRhoMax;
  _srealFunctionMap["RhoMin"] = &Parameters::setRhoMin;
  _srealFunctionMap["WeightEpigraphCutInit"] = &Parameters::setWeightEpigraphCutInit;
  _srealFunctionMap["RadiusBall"] = &Parameters::setRadiusBall;
  _srealFunctionMap["EpsilonReal"] = &Parameters::setEpsilonReal;
  _srealFunctionMap["EpsilonTol"] = &Parameters::setEpsilonTol;
  _srealFunctionMap["WeightEpigraphCutInc"] = &Parameters::setWeightEpigraphCutInc;

  _grealFunctionMap["ObjectiveLB"] = &Parameters::getObjLB;
  _grealFunctionMap["ObjectiveUB"] = &Parameters::getObjUB;
  _grealFunctionMap["Delta"] = &Parameters::getDelta;
  _grealFunctionMap["Eta"] = &Parameters::getEta;
  _grealFunctionMap["Gamma"] = &Parameters::getGamma;
  _grealFunctionMap["Tolerance"] = &Parameters::getTolerance;
  _grealFunctionMap["Rho"] = &Parameters::getRho;
  _grealFunctionMap["RhoMax"] = &Parameters::getRhoMax;
  _grealFunctionMap["RhoMin"] = &Parameters::getRhoMin;
  _grealFunctionMap["WeightEpigraphCutInit"] = &Parameters::getWeightEpigraphCutInit;
  _grealFunctionMap["RadiusBall"] = &Parameters::getRadiusBall;
  _grealFunctionMap["EpsilonReal"] = &Parameters::getEpsilonReal;
  _grealFunctionMap["EpsilonTol"] = &Parameters::getEpsilonTol;
  _grealFunctionMap["WeightEpigraphCutInc"] = &Parameters::getWeightEpigraphCutInc;

  // The string values functions
  _sstringFunctionMap["ProblemName"] = &Parameters::setProblemName;
  _sstringFunctionMap["ProblemAcronym"] = &Parameters::setProblemAcronym;
  _sstringFunctionMap["MethodName"] = &Parameters::setMethodName;
  _sstringFunctionMap["OptimizationType"] = &Parameters::setOptimizationType;
  _sstringFunctionMap["LPSolverName"] = &Parameters::setLPSolverName;

  _gstringFunctionMap["ProblemName"] = &Parameters::getProblemName;
  _gstringFunctionMap["ProblemAcronym"] = &Parameters::getProblemAcronym;
  _gstringFunctionMap["MethodName"] = &Parameters::getMethodName;
  _gstringFunctionMap["OptimizationType"] = &Parameters::getOptimizationType;
  _gstringFunctionMap["LPSolverName"] = &Parameters::getLPSolverName;
} 

bool
Parameters::addParameter(const string &paramName, const string &paramType, const string &paramValue)
{
  if (paramType == "I" || paramType == "Int") {
    std::istringstream instream(paramValue);
    int value;
    instream >> value;
    bool result = setIntParameter(paramName.c_str(), value);
    return result;
  } 
  
  if (paramType == "D" || paramType == "Double" || paramType == "F" || paramType == "Float") {
    std::istringstream instream(paramValue);
    double value;
    instream >> value;
    bool result = setRealParameter(paramName.c_str(), value);
    return result;
  }

  if (paramType == "S" || paramType == "String") {
    bool result = setStringParameter(paramName.c_str(), paramValue.c_str());
    return result;
  }
  AccpmWarning("Unknown parameter type.");
  return false;
}

Parameters::Parameters()
{
  init();
}

Parameters::Parameters(const char *fileName)
{
  std::ifstream fin(fileName);
  if (!fin.is_open()) {
    AccpmError("Cannot open file: ");
    std::cout << fileName << std::endl;
    exit(1);
  }
  
  init();
  string line;
  string paramName;
  string paramType;
  string paramValue;
  
  std::istringstream instream;
  int lineCount = 0;
  while (getline(fin, line)) {
    ++lineCount;
    instream.clear();     // Reset from possible previous errors.
    instream.str(line); 
    if ((instream >> paramName >> paramType >> paramValue >> std::ws) && instream.eof()) {
      if (paramName.find("#") == 0) { // Ignore lines starting with "#"
	//cout << "Ignoring line" << lineCount << endl;
      } else {
	//cout << "Read" << paramName << "|" << paramType << "|" << paramValue << endl;
	if (!addParameter(paramName, paramType, paramValue)) {
	  AccpmWarning("Ignoring line: "); 
	  cout << lineCount << endl;
	}
      }
    } else {
      if (paramName.find("#") == 0 || paramName.find("\n")) {
	//cout << "Ignoring line" << lineCount << endl;
	continue;
      } 
      cout << "Bad Input" << endl;
    }
  }
}

Parameters::~Parameters()
{
  delete _startingPoint;
  delete _varLB;
  delete _varUB;
  delete _pi;
  delete _b;
  delete _centerBall;
  delete _avIndex;
  removeEqualityConstraints();

}

bool
Parameters::setIntParameter(const char *name, int value)
{
  if (_sintFunctionMap.find(name) != _sintFunctionMap.end()) {
    return ((this->*_sintFunctionMap[name])(value));
  } else {
    cerr << "Error in Parameters::setIntParameter : Parameter " << name << " does not exist" << endl;
    return false;
  } 

}

int
Parameters::getIntParameter(const char *name) const
{
  GetIntFunctorMap::const_iterator iter = _gintFunctionMap.find(name);
  if (iter != _gintFunctionMap.end()) {
    return ((this->*((*iter).second))());
  } else {
    cerr << "Error in Parameters::getIntParameter : Parameter " << name << " does not exist" << endl;
    return 0;
  } 
}

bool
Parameters::setRealParameter(const char *name, Real value)
{
  if (_srealFunctionMap.find(name) != _srealFunctionMap.end()) {
    return ((this->*_srealFunctionMap[name])(value));
  } else {
    cerr << "Error in Parameters::setRealParameter : Parameter " << name << " does not exist" << endl;
    return false;
  } 

}

Real
Parameters::getRealParameter(const char *name) const
{
  GetRealFunctorMap::const_iterator iter = _grealFunctionMap.find(name);
  if (iter != _grealFunctionMap.end()) {
    return ((this->*((*iter).second))());
  } else {
    cerr << "Error in Parameters::getRealParameter : Parameter " << name << " does not exist" << endl;
    return 0;
  } 
}
  
bool
Parameters::setStringParameter(const char *name, const string &value)
{
  if (_sstringFunctionMap.find(name) != _sstringFunctionMap.end()) {
    return ((this->*_sstringFunctionMap[name])(value));
  } else {
    cerr << "Error in Parameters::setStringParameter : Parameter " << name << " does not exist" << endl;
    return false;
  } 

}

const string
Parameters::getStringParameter(const char *name) const
{
  GetStringFunctorMap::const_iterator iter = _gstringFunctionMap.find(name);
  if (iter != _gstringFunctionMap.end()) {
    return ((this->*((*iter).second))());
  } else {
    cerr << "Error in Parameters::getStringParameter : Parameter " << name << " does not exist" << endl;
    return string();
  } 
}

bool 
Parameters::setOptimizationType(const string &value)
{
  if (value == "Max") {
    _optimizationType = value;
    _type = OPT_MAX;
  } else {
    if (value == "Min") {
      _optimizationType = value;
      _type = OPT_MIN;
    } else {
      AccpmWarning("Parameters::setOptimizationType: Only supported values are 'Max' and 'Min'");
      AccpmWarning("Setting value to 'Min'");
    }
  }
  return true;
}
    
const string&
Parameters::getOptimizationType(void) const
{
  return _optimizationType;
}

bool
Parameters::setMethodName(const string &value)
{
  if (value == "ProximalACDual") {
    _methodName = value;
    return true;
  } else {
    AccpmError("Currently only supported method is ProximalACDual");
    _methodName = "ProximalACDual";
    return false;
  }
}

bool
Parameters::setLPSolverName(const string &value)
{
  bool status = true;
  if (value == "GLPK") {
    _lpSolverName = "GLPK";
  }
  else {
    if (value == "CLP") {
      _lpSolverName = "CLP";
    } else {
      AccpmError("Only GLPK/CLP LP solvers supported. Using GLPK solver.");
      _lpSolverName = "GLPK";
      status = false;
    }
  }

  return status;
}

bool 
Parameters::setStartingPoint(const StdRealVector &v)
{
  if (_numVariables != v.size()) {
    AccpmError("In Parameters::setStartingPoint():  : size of vector should be equal to number of variables");
    return false;
 } else {
    if (_startingPoint) {
      delete _startingPoint;
    } 
    _startingPoint = new AccpmVector(v.size());
    _startingPoint->copy(v);
    if (_ball) {
      setCenterBall(v);
    }
    return true;
  }
}

bool 
Parameters::setVariableLB(const StdRealVector &v)
{
  if (_numVariables != v.size()) {
    AccpmError("In Parameters::setVariableLB(): size of vector should be equal to number of variables");
    return false;
  } else {
    if (_varLB) {
      delete _varLB;
    } 
    _varLB = new AccpmVector(v.size());
    _varLB->copy(v);
    return true;
  }
}

bool 
Parameters::setVariableUB(const StdRealVector &v)
{
  if (_numVariables != v.size()) {
    AccpmError("In Parameters::setVariableUB(): size of vector should be equal to number of variables"); 
    return false;
  } else {
    if (_varUB) {
      delete _varUB;
    } 
    _varUB = new AccpmVector(v.size());
    _varUB->copy(v);
    return true;
  }
}

bool 
Parameters::setB(const StdRealVector &v)
{
  if (_numVariables != v.size()) {
    AccpmError("In Parameters::setB(): size of vector should be equal to number of variables"); 
    return false;
  } else {
    if (_b) {
      delete _b;
    } 
    _b = new AccpmVector(v.size());
    _b->copy(v);

    return true;
  }
}

bool 
Parameters::setB(const AccpmVector &v)
{
  if ((int)_numVariables != v.size()) {
    AccpmError("In Parameters::setB(): size of vector should be equal to number of variables"); 
    return false;
  } else {
    if (_b) {
      delete _b;
    } 
    _b = new AccpmVector(v.size());
    _b->RealVector::copy(v);
 
    return true;
  }
}

bool 
Parameters::setPi(const StdRealVector &v)
{
  if (_numSubProblems != v.size()) {
    AccpmError("In Parameters::setPi(): size of vector should be equal to number of subproblems "); 
    return false;
  } else {
    if (_pi) {
      delete _pi;
    } 
    _pi = new AccpmVector(v.size());
    _pi->copy(v);
    return true;
  }
}

bool 
Parameters::setCenterBall(const StdRealVector &v)
{
  if (_numVariables != v.size()) {
    AccpmError("In Parameters::setCenterBall(): size of vector should be equal to number of variables"); 
    return false;
  } else {
    if (_centerBall) {
      delete _centerBall;
    } 
    _centerBall = new AccpmVector(v.size());
    _centerBall->copy(v);
    return true;
  }
}

bool
Parameters::setBall(int n)
{
  n == 1 ? _ball = true : _ball = false; 
  if (_ball && _centerBall == 0) {
    if (_startingPoint) {
      _centerBall = new AccpmVector(_numVariables);
      _centerBall->RealVector::copy(*_startingPoint);
    } else {
      setCenterBall(StdRealVector(_numVariables, 0));
    }
  }
  return true; 
}

bool
Parameters::setDiagHessian(int n)
{
  _diagHessian = true;
  if (n != 1) {
    AccpmWarning("In Parameters::setDiagHessian(): currently only Diagonal Hessian is supported. Keeping parameter DiagHessian to 1");
    return false;
  } 
  return true;
}

bool
Parameters::setOracle(Oracle *oracle)
{
  _oracle = oracle;
  return true;
}

const Oracle*
Parameters::getOracle(void) const
{
  return _oracle;
}

/**
 * Set the vector of indices of variables which are active.
 */
bool 
Parameters::setActiveVariableIndex(const StdIntSet &v)
{
  if (_avIndex) {
    delete _avIndex;
    _avIndex = 0;
  }
  if (!v.empty()) {
    _avIndex = new StdIntSet(v);
  } else {
    AccpmError("Deactivating all variables");
  }

  return true;
}

/**
 * Get the vector of indices of variables which are active.
 */
const StdIntSet *
Parameters::getActiveVariableIndex() const
{
  return _avIndex;
}

/**
 * Is the Active Set strategy supported.
 */
bool 
Parameters::hasVariableDimension() const
{
  if (_avIndex == 0) {
    return false;
  } else {
    return true;
    // return _avIndex->size() == _numVariables ? false : true;
  }
}

bool 
Parameters::setProximal(int n) 
{ 
  n == 1 ? _proximal = true : _proximal = false; 
  _rho = n;
  if (n == 0) {
    AccpmWarning("Proximal term is being turned off. Are you sure?\n");
  }
  return true; 
}

int 
Parameters::getProximal() const 
{ 
  return _proximal; 
}

void 
Parameters::addEqualityConstraints(const AccpmGenMatrix &constraints, const AccpmVector &rhs)
{
  assert(constraints.size(0) == (int)_numVariables);
  assert(constraints.size(1) == rhs.size());
  int numConstraints = constraints.size(1);
  int previousNumConstraints = 0;
  if (_D == 0) {
    _D = new AccpmGenMatrix(constraints);
    _d = new AccpmVector;
  } else {
    previousNumConstraints = _D->size(1);
    numConstraints += previousNumConstraints;
    AccpmGenMatrix tmp(_D->size(0), numConstraints);
    for (int i = 0; i < previousNumConstraints; ++i) {
      tmp.assignColumn(i, _D->getColumn(i));
    }
    _D->ref(tmp);
  }
  for (int i = 0; i < constraints.size(1); ++i) {
    _D->assignColumn(previousNumConstraints + i, constraints.getColumn(i));
  }
  _d->append(rhs);
  std::cout << "Added constraints:\n" << *_D << "\n" << *_d << std::endl;
}

void
Parameters::removeEqualityConstraints()
{
  delete _D;
  _D = 0;
  delete _d;
  _d = 0;
}

bool
Parameters::hasEqualityConstraints() const
{
  return _D != 0;
}

void
Parameters::getEqualityConstraints(AccpmGenMatrix *&constraints, AccpmVector *&rhs) const
{
  constraints = _D;
  rhs = _d;
}

void 
Parameters::output(std::ostream &os) const
{
  
  os << "\n-------------------------------------------------------------------------" << endl;
  os << "Parameters: " << endl;
  os << "Problem Specific Parameters for: "  << _problemName << "(" << _problemAcronym << ")" << endl;
  os << "Optimization Type:\t\t" << getOptimizationType() << endl;
  os << "Number of Variables:\t\t" << _numVariables << endl;
  os << "Number of Subproblems:\t\t" << _numSubProblems << endl;
  
  os << "Objective Lower Bound:\t\t" << _objLB << endl;
  os << "Objective Upper Bound:\t\t"  << _objUB << endl;
  os << endl;
  os << "\nMethod Parameters:" << endl;
  os << "Method Name:\t\t" << _methodName << endl;
  os << "Delta:\t\t\t" << _delta << endl;
  os << "Eta:\t\t\t" <<  _eta << endl;;
  os << "Gamma:\t\t\t" << _gamma << endl;
  os << "Tolerance:\t\t" << _tolerance << endl;
  os << "Max Outer Iterations: " << _maxOuterIterations 
     << " Max Inner Iterations:" << _maxInnerIterations << endl;
  os << "Proximal:\t\t" << _proximal << endl;
  os << endl;
  os << "Filter:\t\t" << _filter << endl;
  os << "Convexity Check: " << _convexityCheck 
     << " ConvexityFix: " << _convexityFix << endl;
  os << endl;
  os << "Rho:\t\t" << _rho << endl;
  os << "DynamicRho:\t\t" << _dynamicRho << endl;
  os << "Rho Max: " << _rhoMax 
     << " RhoMin: " << _rhoMin << endl;
  os << endl;
  os << "WeightEpigraphCutInit:\t" << _weightEpigraphCutInit << endl;
  os << "WeightEpigraphCutInc:\t" << _weightEpigraphCutInc << endl;
  os << "Ball:\t" << _ball << " Radius Ball:\t" << _radiusBall << endl;

  os << "-------------------------------------------------------------------------" << endl;
  if (_verbosity > 1) {
    if (_startingPoint) {
      os << "Starting Point:" << "\n" << *_startingPoint; 
    }
    if (_varLB) {
      os << "Variable LB:" << "\n" << *_varLB;
    }
    if (_varUB) {
      os << "Variable UB:" << "\n" << *_varUB;
    }
    if (_b) {
      os << "Linear Component b:" << "\n" << *_b;
    }
    if (_pi) {
      os << "Weight Vector pi:" << "\n" << *_pi;
    }
    if (_centerBall) {
      os << "Center Ball:" << "\n" << *_centerBall;
    }
  }
}

std::ostream& 
operator<<(std::ostream &os, const Parameters &P)
{
  P.output(os);
  return os;
}

}
