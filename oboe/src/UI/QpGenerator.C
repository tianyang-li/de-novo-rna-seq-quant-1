// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "QpGenerator.h"
#include "Oracle.h"
#include "Parameters.h"
#include "Manager.h"
#include "Method.h"
#include "AccpmBlasInterface.h"
#include <fstream>
#include <iomanip>

using std::endl;
namespace Accpm {

QpGenerator::QpGenerator() : _manager(0), _method(0)
{
  _timer.start();
  _oracleTime = 0;
}

QpGenerator::~QpGenerator()
{
  _param = 0;
  delete _manager;
  delete _method;
}

void
QpGenerator::init(Parameters *param, Oracle *oracle)
{
  std::cout << "************************************************************************" << std::endl;
  std::cout << "OBOE Version 1.0" <<
    "\n(C) 2005 Université de Genève. See Copyright Notification in OBOE" << std::endl;
  std::cout << "************************************************************************\n" << std::endl;
  _param = param;
  _param->setOracle(oracle);
  if (_param->getStringParameter("MethodName") == "ProximalACDual") {
    _method = new DualMethod(_param);
  }
  _manager = new Manager(_param);
}

int
QpGenerator::run()
{
  if (_manager->getNumOuterIteration() < _param->getIntParameter("MaxOuterIterations")) {
    const AccpmVector &y = _manager->getCurrentY();
    int n = y.size();
    assert(n ==  _param->getIntParameter("NumVariables"));
    int p = _param->getIntParameter("NumSubProblems");
    AccpmGenMatrix subGrad(n, p);
    AccpmVector val(p);
    AccpmGenMatrix subProblemIndex(p,1);
    if (_param->getIntParameter("Verbosity") > 0) {
      printHeader();
    }

    if (_param->getOracle()) {
      OracleFunction *f1 = _param->getOracle()->getF1();
      assert(f1);
      Timer t;
      t.start();
      if (!f1->eval(y, val, subGrad, &subProblemIndex)) {
	t.stop();
	_oracleTime  += t.getCpuTime();
	if (_param->getOracle()->computesBound()) {
	  double lb = _param->getOptType() * _param->getOracle()->getObjectiveFunctionBound();
	  lb = std::max(lb, _manager->getObjLB());
	  _manager->updateLB(lb);
	  
	}
	_manager->update1(val, subGrad, subProblemIndex);
	
	OracleFunction *f2 = _param->getOracle()->getF2();
	AccpmVector val2(1);
	val2 = 0;
	if (_manager->isCurrentPointFeasible() && f2) {
	  AccpmGenMatrix subGrad2(n, 1);
	  AccpmGenMatrix hessian;
	  if (_param->getIntParameter("DiagHessian")) {
	    hessian.resize(n, 1);
	  } else {
	    hessian.resize(n, n);
	  }
	  t.start();
	  f2->eval(y, val2, subGrad2, &hessian);
	  t.stop();
	  _oracleTime  += t.getCpuTime();
	  _manager->update2(val2(0), subGrad2, hessian);
	}
	if (_param->getIntParameter("Verbosity") > 0) {
	  printIteration(val, val2, subGrad, subProblemIndex);
	}
	_method->run(*_manager);	
        _manager->updatePreviousPhase2();
	if (_manager->getExitCode() == RELATIVE_GAP_REACHED) {
	  std::cout << "Relative Gap reached" << std::endl;
	  return 1;
	}
      } else {
	t.stop();
	_oracleTime  += t.getCpuTime();
	std::cout << "User Oracle Stop" << std::endl;
	_manager->setExitCode(USER_STOP);
	return 3;
      }
    }
    return 0;
  }
  _manager->setExitCode(MAX_OUTER_ITERATIONS);
  std::cout << "Maximum number of outer iterations exceeded" << std::endl;
  return 1;
}

void
QpGenerator::printHeader(std::ostream &os) const
{
  if (_manager->getNumOuterIteration() == 0) {
    if (_param->getIntParameter("Verbosity") == 1) {
      os << "Iteration\tObjective" << std::endl;
    } else if (_param->getIntParameter("Verbosity") == 2) {
      os << "Outer\tInner\tFeasible\t\tObjective\t\tUBound\t\tLBound\tRelativeGap" << std::endl;	
    }
  }
}

void
QpGenerator::printIteration(const AccpmVector &val, const AccpmVector &val2, const AccpmGenMatrix &subGrad,
			    const AccpmGenMatrix &subProblemIndex, std::ostream &os) const
{
  const AccpmVector &y = _manager->getCurrentY();
  double fVal = val2(0); 
  if (_manager->isCurrentPointFeasible()) { 
    fVal += AccpmLADotProd(val, *_param->getPi()); 
    if (_param->getB()) {
      fVal +=  AccpmLADotProd(y, *_param->getB());
	  }
  } else {
    fVal += val.sum();  /* For now we display the sum of cut values for fesibility cut */ 
  }
  if (_param->getIntParameter("Verbosity") > 2) {
    std::cout << "Iteration " << _manager->getNumOuterIteration() << ":" << std::endl;
    std::cout << "y:\n " << y << std::endl;
    std::cout << "Objective:" << fVal << "\n" << std::endl;
    if (_param->getIntParameter("Verbosity") > 3) {
      std::cout << "Subgradient:\n" << subGrad << std::endl;
      std::cout << "SubProblemIndex:\n" << subProblemIndex << std::endl;
    }
  } else {
    if (_param->getIntParameter("Verbosity") == 1) {
      std::cout << std::setw(4) << _manager->getNumOuterIteration() 
		<< std::setw(18) << fVal << "\n" << std::endl;
    } else if (_param->getIntParameter("Verbosity") == 2) {
      std::cout << std::setw(4) << _manager->getNumOuterIteration() 
		<< std::setw(8) << _manager->getNumInnerIterations() 
		<< std::setw(8) << _manager->isCurrentPointFeasible()
		<< " " << std::setw(16) << fVal 
		<< " " << std::setw(16) << getObjUB() 
		<< " " << std::setw(16) << getObjLB() 
		<< "\t" << _manager->getRelativeGap() << std::endl;  
    }
  }
}

void 
QpGenerator::terminate()
{
}

void
QpGenerator::output(std::ostream &os)
{
  _timer.stop();
  if (_manager) {
    os << "\n-------------------------------------------------------------------------" << endl;
    os << "OBOE using ACCPM - Vital Statistics" << endl;
    _manager->output(os);
    os << "Number of Outer Iterations:\t" << _manager->getNumOuterIteration() << endl;
    os << "Inner / Outer Iterations:\t" << 
      _manager->getNumInnerIterations()*1.0/_manager->getNumOuterIteration() << endl;
    double totalTime = _timer.getCpuTime();
    double oraclePercent = 0;
    if (totalTime > 0) {
      oraclePercent = (_oracleTime * 100) / totalTime;
    }
    os << "Total Elapsed Time:\t\t" << _timer.getRealTime() << " sec" << endl;
    os << "Total Cpu Time:\t\t\t" << totalTime << " sec" << endl;
    os << "Oracle Time:\t\t\t" << _oracleTime << " sec (" 
       << oraclePercent << " %)" 
       << endl;

    if (_param->getIntParameter("Verbosity") > 5) {
      os << "\nSolution y:\n" << _manager->getBestY() << std::endl;
      os << "Input Parameters:" << *_param << std::endl;
    }
    os << "\n-------------------------------------------------------------------------" << endl;
  }
}

const AccpmVector *
QpGenerator::getQueryPoint() const
{
  return &_manager->getCurrentY();
}

double 
QpGenerator::getOptimalObj() const
{
  return _param->getOptType() * _manager->getObjUB();
}

double 
QpGenerator::getObjUB() const
{
  if (_param->getOptType() == OPT_MAX) {
    return _param->getOptType() * _manager->getObjLB();
  }
  return _manager->getObjUB();
}

double 
QpGenerator::getObjLB() const
{
  if (_param->getOptType() == OPT_MAX) {
    return _param->getOptType() * _manager->getObjUB();
  }
  return _manager->getObjLB();
}

double 
QpGenerator::getRelativeGap() const
{
  return _manager->getRelativeGap();
}

int 
QpGenerator::getNumCuts() const
{
  return _manager->getNumCuts();
}

int
QpGenerator::getActiveCuts(AccpmGenMatrix &cuts) const
{
  return _manager->getActiveCuts(cuts);
}

int
QpGenerator::getActiveCuts(const AccpmGenMatrix *&cuts) const
{
  return _manager->getActiveCuts(cuts);
}

const AccpmVector&
QpGenerator::getRhs() const
{
  return _manager->getRhsCoef();
}

const AccpmVector& 
QpGenerator::getCurrentZ() const
{
  return _manager->getCurrentZ();
}

const AccpmVector& 
QpGenerator::getCutType() const
{
  return _manager->getSubProblemIndex();
}

const AccpmVector *
QpGenerator::getCurrentX() const
{
  const AccpmVector *x = &(_manager->getCurrentX());
  if (x->size() > 0) {
    return x;
  }
  return 0;
}

ExitCode 
QpGenerator::getExitCode() const
{
  return _manager->getExitCode();
}

#ifdef SERIALIZATION
template<class Archive> 
void
QpGenerator::serialize(Archive &ar, const unsigned int file_version)
{
  std::cout << "Serializing QpGenerator..." << std::endl;
  ar & *_manager;
}
#endif

void
QpGenerator::save(const char *fileName) const 
{
#ifdef SERIALIZATION
  std::ofstream ofs(fileName);
  boost::archive::text_oarchive oa(ofs);
  oa << *this;
#else
  std::cout << "Error: Serialization libraries are not present. Cannot save/load." 
	    << std::endl;
#endif
}

void
QpGenerator::load(const char *fileName) 
{
#ifdef SERIALIZATION
  // open the archive
  std::ifstream ifs(fileName);
  boost::archive::text_iarchive ia(ifs);
  
  // restore the schedule from the archive
  ia >> *this;
#else
  std::cout << "Error: Serialization libraries are not present. Cannot save/load." 
	    << std::endl;
#endif
}

}
