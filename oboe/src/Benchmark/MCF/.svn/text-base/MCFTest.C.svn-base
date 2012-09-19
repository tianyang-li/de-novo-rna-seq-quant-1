// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// The OBOE team
//

#undef LA_BOUNDS_CHECK

#include "QpGenerator.h"
#include "Parameters.h"
#include "MCFOracle.h"
#include "MCFOracleFunction.h"
#include "MCFSmoothOracleFunction.h"
#include "AccpmBlasInterface.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>


using std::string;

static int
getLineCount(const char *fileName)
{
  int lc = 0;
  std::ifstream fin(fileName);
  string line;
  if (fin.is_open()) {
    while (getline(fin, line)) {
      ++lc;
    }
  } else {
    AccpmError("Cannot open file: ");
    std::cout << fileName << std::endl;
    exit(1);
  }
  fin.close();
  return lc;
}

int 
main(int argc, char *argv[])
{
  char *paramFile = "param.txt";
  char *networkFile = "Cpl30.txt";
  char *demandFile = "Dpl30.txt";
  char type = 'l';
  if (argc > 1) {
    if (argc != 4 && argc != 5) {
      std::cout << "Error starting " << argv[0] << std::endl;
      std::cout << "Usage: " << argv[0] << " <Parameter File> <Network Data File>"
		<< " <Demand Data File> [<Function Type: Linear/Kleinrock>]" 
		<< std::endl;
      exit(0);
    }
    paramFile = argv[1];
    networkFile = argv[2];
    demandFile = argv[3];
    if (argc == 5) {
      type = *argv[4];
    }
  } 
  int n = getLineCount(networkFile)/2;
  int numK = getLineCount(demandFile)/2;
  
  Accpm::Parameters param(paramFile);
  if (n != param.getIntParameter("NumVariables")) {
    param.setIntParameter("NumVariables", n);
  }
  
  MCFOracleFunction f1(networkFile, demandFile, n, numK, type);
  OracleFunction *f2;
  if (type == 'l') {
      f2 = 0;
      param.setB(f1.getCapacity());
  } else {
    f2 = new MCFSmoothOracleFunction(f1.getMCFData());
  }
  MCFOracle oracle(&f1, f2, f1.getMCFData());
  
  int p = f1.getMCFData()->getNumSources();
  if (param.getIntParameter("NumSubProblems") > 1) {
    if (p !=  param.getIntParameter("NumSubProblems")) {
      std::cout << "Function is being disaggregated. Number of subproblems modified to : " << p << std::endl;
      param.setIntParameter("NumSubProblems", p);
    }
    f1.setAggregationLevel(2);
  } else {
      p = 1;
  }
  
  vector<double> start(n, 1e-3);
  vector<double> varLB(n, 0);
  if (tolower(type) == 'k') {
    AccpmVector c3(n);
    c3 = 1;
    c3.rdivide(f1.getCapacity());
    for (int i = 0; i < n; ++i) {
      start[i] += c3(i);
      varLB[i] += c3(i);
    }
  }
  param.setStartingPoint(start);
  param.setVariableLB(varLB);
  
  vector<double> centerBall(n, 3);
  param.setCenterBall(centerBall);
    
    
  QpGenerator qpGen;
  qpGen.init(&param, &oracle);
  
  std::cout << "Solving MCF for " << f1.getNumNodes() << " nodes, " 
	    << n << " arcs and " << numK << " commodities" << std::endl;
  
  AccpmVector Ax(n); 
  const AccpmGenMatrix *cuts;
  const AccpmVector *x;
  //while (DBL_GT(qpGen.getRelativeGap(), param.getRealParameter("Tolerance"))) {
   while(!qpGen.run()) {
    if (qpGen.getExitCode() == CHOLESKY_FAILURE) {
	break;
    }
    qpGen.getActiveCuts(cuts);
    x = qpGen.getCurrentX();
      
    if (x) {
      AccpmVector xScaled(*x);
      double factor = (p*1.0)/x->sum();
   //   std::cout << "x->sum()" << x->sum() << std::endl;
      AccpmLAScale(factor, xScaled);
      double sum = xScaled.sum();
      if (!DBL_CMP(sum, p, 1e-12)) {
	std::cerr << std::setprecision(16) 
		  << "xScaled not summing to 1: " 
		  << std::setprecision(16) << sum << std::endl;
      }
      AccpmLAMatVecMult(*cuts, xScaled, Ax, -1, 0);
      f1.updateAx(Ax, factor);
      oracle.computeLowerBound();
    }
  }
  x = qpGen.getCurrentX();
  if (x) {
    qpGen.getActiveCuts(cuts);
    AccpmVector xScaled(*x);
    double factor = (p*1.0)/x->sum();
    AccpmLAScale(factor, xScaled);
    std::cout << "Sum of x: " << x->sum() << std::endl;
    std::cout << "Sum of xScaled: " << xScaled.sum() << std::endl;
    //  std::cout << "x:\n" << *x << std::endl;
    AccpmLAMatVecMult(*cuts, xScaled, Ax, -1, 0);
    int numSatArcs = 0;
    const AccpmVector &capacity = f1.getCapacity();
      for (int i = 0; i < capacity.size(); ++i) {
	if (!DBL_LT(Ax(i), 0.99*capacity(i))) {
	  ++numSatArcs;
	}
      }
      std::cout << "Number of Saturated Arcs: " << numSatArcs  << "(" 
		<< (numSatArcs * 100.0)/n << " %)" << std::endl;
  }
  qpGen.output(std::cout);
  
  qpGen.terminate();
  //getchar();
}

