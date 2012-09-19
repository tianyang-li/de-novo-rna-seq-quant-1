// Copyright (c) 2004-2007 University de Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "Oracle.h"
#include "QpGenerator.h"
#include "Parameters.h"
#include "AccpmDynMatrix.h"
#include "AccpmBlasInterface.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace Accpm;
using namespace std;
using std::string;

static void
printUsage(const char *progName)
{  
  std::cout << "Usage: " << progName << " <Parameter File> <Input Data File>"
	    << " <Number of Attributes(including the binary label> <Number of Points> <nu>" 
	    << std::endl;
}

/**
 * Oracle for using the homogeneous version of ACCPM (h-ACCPM).
 * This oracle is used to solve the datamining problem which
 * attempts to find a linear separation of points divided into
 * 2 classes.
 *
 * The input Data file is a matrix giving the values for each attribute.
 * Each line of this file is interpreted as an observation(point) with 
 * the values of attributes, the last of which is a flag 0/1 giving the 
 * class of that point.
 */
class hAccpmOracleFunction : public OracleFunction {

private:
  int _numAttributes;
  int _numPoints;
  double _nu;
  
  AccpmGenMatrix _B0;
  AccpmGenMatrix _B1;
  
  void readMatrix(const string &fileName) {
    //std::cout << "Reading File:" << fileName << std::endl;
    std::cout << "Number of Attributes:" << _numAttributes << std::endl;
    std::ifstream fin(fileName.c_str());
    if (fin.is_open()) {
      string line;
      std::istringstream instream;
      double val;
      AccpmDynMatrix B0;
      AccpmDynMatrix B1;
      for (int i = 0; i < _numPoints; ++i) {
	if (getline(fin, line)) {
	  AccpmVector v;
	  instream.clear();
	  instream.str(line); 
	  for (int j = 0; j < _numAttributes; ++j) { 
	  if ((instream >> val >> std::ws)) {
	    v.append(val); 
	  } else {
	    std::cerr << "Error reading matrix from file: " << fileName 
		      << " at line:" << i+1 << " column: " << j + 1 << std::endl;
	  }
	  }
	  if (v(_numAttributes - 1) == 0) {
	    v(_numAttributes - 1) = 1; // putting last row of ones in both B0 and B1 for easier matrix operations
	    B0.addColumn(v);
	  } else if (v(_numAttributes - 1) == 1) {
	    B1.addColumn(v);
	}
	} 
      }
      _B0 = B0.getM();
      _B1 = B1.getM();
      if (_numPoints != (_B1.size(1) + _B0.size(1))) {
	std::cout << "Mismatch in number of points: only " << _B1.size(1) + _B1.size(1)
		  << " of " << _numPoints << " points read" << std::endl;
	std::cout << "The others possibly have missing values" << std::endl;
      _numPoints = _B1.size(1) + _B1.size(1);
      }
      
      std::cout << "Classifcation of " << _numPoints << " points: [" << _B1.size(1) << " , " << _B0.size(1) 
		<< "]" << std::endl;
      // std::cout << "B0:\n" << _B0 << std::endl;
      //std::cout << "B1:\n" << _B1 << std::endl;
    } else {
      std::cout << "\nError opening file: " << fileName << std::endl;
      printUsage("oboeHACCPM");
      exit(1);
    }
  } 
  
public:
  int _numMC0;
  int _numMC1;
  AccpmVector _y;
  hAccpmOracleFunction() : OracleFunction() {}
  hAccpmOracleFunction(const char *fileName, int n, int m, double nu) : OracleFunction(),
									_numAttributes(n), _numPoints(m), _nu(nu) 
  {
    _numMC1 = 0;
    _numMC0 = 0;
    readMatrix(fileName);
  }
  /**
   * h-Accpm Oracle:
   * We need to solve the following optimization problem
   *  Min   f(x) = 1/|S1|sum(max(-w^Tau+gamma+mu,0))+1/|S2|sum(max(w^Tau-gamma+mu,0))
   */
  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {

    _y = y;
    int n = y.size();
    AccpmVector x = y(LaIndex(0, n-2)); 
    AccpmVector subGradx(x.size());
    AccpmLAScale(1.0/y(n-1), x);
    subGradx = 0;
    int n1 = _B1.size(1);
    AccpmVector violation(n1);
    AccpmVector set(n1);
    _numMC1 = 0;
    _numMC0 = 0;
    double fx = 0;
    if  (n1) {
      violation = _nu;
      AccpmLAMatTransVecMult(_B1, x, violation, -1, 1);
      //std::cout << "violation1\n" << violation << std::endl;
      set = 0;
      for (int i = 0; i < n1; ++i) {
	if (violation(i) > 0) {
	  set(i) = 1;
	  _numMC1++;
	}
      }
      
      AccpmLAMatVecMult(_B1, set, subGradx, -1.0/n1, 0);
      fx += (AccpmLADotProd(violation, set) * 1.0) / n1;
    }
    int n0 = _B0.size(1);
    if (n0) {
      violation = AccpmVector(n0);
      violation = _nu;
      
      AccpmLAMatTransVecMult(_B0, x, violation, 1, 1);
      //std::cout << "violation0\n" << violation << std::endl;
      set = AccpmVector(n0);
      set = 0;
      for (int i = 0; i < n0; ++i) {
	if (violation(i) > 0) {
	  set(i) = 1;
	  _numMC0++;
	}
      }
      AccpmLAMatVecMult(_B0, set, subGradx, 1.0/n0, 1);
      fx += (AccpmLADotProd(violation, set) * 1.0) / n0;
    }
    std::cout << "numMC1: " << _numMC1 << " numMC0: " << _numMC0 << std::endl;

    functionValue = 0;

    double lastComp = -AccpmLADotProd(x, subGradx);
    subGradx.append(lastComp);
    memcpy(subGradients.addr(), subGradx.addr(), sizeof(double)*subGradx.size());

    *info = 0;

    if (_numMC1 + _numMC0 == 0) {
      return 1;
    } 
    return 0;
  }
};

int 
main(int argc, char *argv[])
{
  int numAttributes = 31;
  int numPoints = 569;
  char *paramFile = "param.txt";
  char *dataFile = "dtm.dat";
  double nu = 1;
  if (argc > 1) {
    if (argc != 6) {
      std::cout << "\nError starting " << argv[0] << std::endl;
      printUsage(argv[0]);
    
      exit(0);
    }
    paramFile = argv[1];
    dataFile = argv[2];
    numAttributes = atoi(argv[3]);
    numPoints = atoi(argv[4]);
    nu = atof(argv[5]);
  } 

  Accpm::Parameters param(paramFile);
  int n = numAttributes+1;
  assert(n  == param.getIntParameter("NumVariables"));
  vector<double> start(n, 1);
  param.setStartingPoint(start);
  vector<double> varLB(n, ACCPM_MINUS_INF);
  varLB[n-1] = 1e-6;
  param.setVariableLB(varLB);

  hAccpmOracleFunction f1(dataFile, numAttributes, numPoints, nu);
  Oracle oracle(&f1);
  
  QpGenerator qpGen;
  qpGen.init(&param, &oracle);
  std::cout << "\nLinear Separation Problem with nu:" << nu << std::endl << std::endl;
  
  int bestMC1 = numPoints;
  int bestMC0 = numPoints;
  AccpmVector bestY; 
  
  while (!qpGen.run()) {
    std::cout << "Misclassified Points: [" << f1._numMC1 << " , " << f1._numMC0 << "]" << std::endl;
    if (bestMC1 + bestMC0 > f1._numMC1 + f1._numMC0) {
      bestMC1 = f1._numMC1;
      bestMC0 = f1._numMC0;
      bestY = f1._y;
    }
  }
  if (bestMC1 + bestMC0 > f1._numMC1 + f1._numMC0) {
      bestMC1 = f1._numMC1;
      bestMC0 = f1._numMC0;
      bestY = f1._y;
    }
  qpGen.output(cout);
  qpGen.terminate();
  std::cout << "Misclassified Points: [" << bestMC1 << " , " << bestMC0 << "]" << std::endl;
  std::cout << "bestY:\n" << bestY << std::endl;
}

