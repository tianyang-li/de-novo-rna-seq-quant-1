// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
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
#include "AccpmBlasInterface.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace Accpm;
using namespace std;

class PMEDOracleFunction : public OracleFunction {

private:
  int _p;
  AccpmGenMatrix _distanceMat;

public:
  PMEDOracleFunction(int p, int numC, int numF, const char *distFileName) 
    : OracleFunction(), _p(p), _distanceMat(numF, numC) 
  {
    std::ifstream fin(distFileName);
    if (fin.is_open()) {
      string line;
      std::istringstream instream;
      double dist;
      for (int i = 0; i < numC; ++i) {
	if (getline(fin, line)) {
	  instream.clear();
	  instream.str(line); 
	  for (int j = 0; j < numF; ++j) { 
	    if ((instream >> dist >> std::ws)) {
	      _distanceMat(j, i) = dist; //Storing the transpose for easier access
	  } else {
	      std::cerr << "Error reading matrix from file: " << distFileName 
			<< " at line:" << i+1 << " column: " << j + 1 << std::endl;
	    }
	  }
	} else {
	  std::cerr << "Error reading matrix from file: " << distFileName << std::endl;
	}
      }
    } else {
      std::cerr << "Error opening file: " << distFileName << std::endl;
      exit(1);
    }
  }

  
  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {
    int dimOfy = y.size() - 1;
    double totalfVal = 0;
    int *point = new int[dimOfy];
    AccpmVector subGrad(y.size());
    subGrad = 0;
    AccpmVector ny = y(LaIndex(0, dimOfy-1));
    ny.negate();
    for (int i = 0; i < dimOfy; i++) {
      int ind = 0;
      double fVal = -y(dimOfy);
      AccpmVector val = ny;
      AccpmLAAddMult(val, 1, _distanceMat.getColumn(i));
      for (int j = 0; j < dimOfy; ++j) {
	if ((val(j) < 0) || i == j)  {
	  point[ind++]= j;
	  fVal += val(j);
	}
      }
      if (fVal < -1e-8) {
	totalfVal -= fVal;
	for (int j = 0; j < ind; ++j){
	  ++subGrad(point[j]);
	}
	++subGrad(dimOfy);
      }
    }
    delete [] point;

    AccpmVector tmp(y.size());
    tmp = 1;
    tmp(dimOfy) = _p; // Number of facilities
    functionValue = -totalfVal + AccpmLADotProd(tmp, y);
    
    AccpmLAAddMult(tmp, -1, subGrad);
   
    memcpy(subGradients.addr(), tmp.addr(), sizeof(double)*tmp.size());
    
    if (info) {
      *info = 1; // Optimality Cut
    }

    return 0; 
  }
};

int 
main(int argc, char *argv[])
{
  int n = 50;
  char *paramFile = "param.txt";
  char *distFileName = "dist.txt";
  bool check = false;
  if (argc > 1) {
    if (argc != 4 && argc != 5) {
      std::cout << "Error starting " << argv[0] << std::endl;
      std::cout << "Usage: " << argv[0] << " <Parameter File> <Number of variables> <Distance Matrix> [check]" 
	<< std::endl;
      exit(0);
    }
    paramFile = argv[1];
    n = atoi(argv[2]);
    distFileName = argv[3];
    if (argc == 5) {
      (strcmp(argv[4], "check") == 0) ? check = true : check = false;
    } 
  } 
  Accpm::Parameters param(paramFile);
  assert(n == param.getIntParameter("NumVariables"));
  param.setOptimizationType("Max");

  vector<double> start(n, 1);
  param.setStartingPoint(start);
  
  vector<double> varLB(n, 0);
  varLB[n-1] = -10;
  param.setVariableLB(varLB);

  int p = std::max( (n-1) / 10, 2);
  PMEDOracleFunction f1(p, n-1, n-1, distFileName);
  Oracle oracle(&f1);
  
  QpGenerator qpGen;
 
  qpGen.init(&param, &oracle);
  std::cout << "Solving p-Median problem for " << p << " medians" << std::endl;
  
  while (!qpGen.run()) {
  }
  
  qpGen.output(cout);
  
  double obj1 = qpGen.getOptimalObj();
  double relativeGap1 = qpGen.getRelativeGap();
  
  qpGen.terminate();

  if (check) {
    param.setIntParameter("MaxOuterIterations", 17);
    qpGen.init(&param, &oracle);
    while (!qpGen.run()) {
    }
    qpGen.output(cout);
    qpGen.save( "PMED.state1");
    qpGen.terminate();

    QpGenerator qpGen2;
    qpGen2.init(&param, &oracle);
    qpGen2.load( "PMED.state1");
    while (!qpGen2.run()) {
    }
    qpGen2.output(cout);
    
    double obj2 = qpGen2.getOptimalObj();
    double relativeGap2 = qpGen2.getRelativeGap();

    assert(obj1 == obj2);
    assert(relativeGap1 == relativeGap2);

    qpGen2.terminate();
  }
}
