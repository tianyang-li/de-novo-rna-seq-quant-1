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
#include "AccpmDynMatrix.h"
#include "AccpmBlasInterface.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace Accpm;
using namespace std;
using std::string;

class QPOracleFunction : public OracleFunction {

private:
  AccpmGenMatrix _A;
  AccpmGenMatrix _B;
  AccpmGenMatrix _F;
  AccpmGenMatrix _C;
  AccpmGenMatrix _E;
  void readMatrix(const string &fileName, AccpmGenMatrix &matrix) {
    //std::cout << "Reading File:" << fileName << std::endl;
    std::ifstream fin(fileName.c_str());
    if (fin.is_open()) {
      string line;
      std::istringstream instream;
      double val;
      for (int i = 0; i < matrix.size(0); ++i) {
	if (getline(fin, line)) {
	  instream.clear();
	  instream.str(line); 
	  for (int j = 0; j < matrix.size(1); ++j) { 
	    if ((instream >> val >> std::ws)) {
	      matrix(i, j) = val; 
	    } else {
	      std::cerr << "Error reading matrix from file: " << fileName 
			<< " at line:" << i+1 << " column: " << j + 1 << std::endl;
	    }
	  }
	} else {
	  std::cerr << "Error reading matrix from file: " << fileName << std::endl;
	}
      }
    } else {
      std::cerr << "Error opening file: " << fileName << std::endl;
      exit(1);
    }
  }
  
public:
  QPOracleFunction() : OracleFunction() {}
  QPOracleFunction(const char *filePrefix, int n, int m) : OracleFunction(), _A(n,1), _B(n,1), 
							   _F(n,1), _C(n,m), _E(n,m) {
    
    string fileName(filePrefix);
    fileName = fileName + "A" +".dat";
    readMatrix(fileName, _A);
    fileName = string(filePrefix) + "B" + ".dat";
    readMatrix(fileName, _B);
    fileName = string(filePrefix) + "C" + ".dat";
    readMatrix(fileName, _C);
    fileName = string(filePrefix) + "E" + ".dat";
    readMatrix(fileName, _E);
    fileName = string(filePrefix) + "F" + ".dat";
    readMatrix(fileName, _F);
  }
  
  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {
    int n = _C.size(0);
    int m = _C.size(1);
    bool feasible = true;
    AccpmVector u, v;
    double w;
    AccpmVector F;
    AccpmDynMatrix G;
    AccpmVector C;

    for (int j = 0; j < m; ++j) {
      u = _C.getColumn(j);
      v = _E.getColumn(j);
      w = _F(j, 0);
      
      AccpmVector uyv = y;
      uyv.times(u);
      AccpmLAAddMult(uyv, -1, v);
      double constraintOneGap = AccpmLADotProd(uyv, uyv);
      constraintOneGap -= w;
      if (constraintOneGap > 0) {
	feasible = false;
	F.append(constraintOneGap);
	//G = [G , (2 * u).*(u .* y - v)];
	AccpmVector subGrad = u;
	AccpmLAScale(2, subGrad);
	subGrad.times(uyv);
	G.addColumn(subGrad);
	C.append(0);
      }
    }
  
    if (feasible) {
      v = _A.getColumn(0);
      v.times(y);
      AccpmLAAddMult(v, -1, _B.getColumn(0));
      AccpmVector v2 = v;
      v2.times(v);
      F.append(v2);

      //G = [G , diag((2 * a).*(a .* y - b))];
      AccpmVector diag = _A.getColumn(0);
      AccpmLAScale(2, diag);
      diag.times(v);
      subGradients = AccpmGenMatrix(n, n);
      subGradients = 0;
      for (int i = 0; i < n; subGradients(i,i) = diag(i), ++i);

      for(int i = 0; i < n; diag(i) = i+1, ++i);
      C.append(diag);
    
    } else {
      subGradients = G.getM();
    }
    
    functionValue = F;
    if (info) {
      AccpmGenMatrix CM(C.size(), 1);
      CM.assignColumn(0, C);
      *info = CM;
    }
    return 0;
  }
};

int 
main(int argc, char *argv[])
{
  int n = 5;
  int m = 5;
  char *paramFile = "param.txt";
  char *matrixFilePrefix = "QP_";
  if (argc > 1) {
    if (argc != 5) {
      std::cout << "Error starting " << argv[0] << std::endl;
      std::cout << "Usage: " << argv[0] << " <Parameter File> <Input Data File prefix>"
		<< " <Number of variables> <Number of constraints>" 
		<< std::endl;
      exit(0);
    }
    paramFile = argv[1];
    matrixFilePrefix = argv[2];
    n = atoi(argv[3]);
    m = atoi(argv[4]);
  } 
  Accpm::Parameters param(paramFile);
  assert(n == param.getIntParameter("NumVariables"));
  vector<double> start(n, 0);
  param.setStartingPoint(start);
  
  QPOracleFunction f1(matrixFilePrefix, n, m);
  Oracle oracle(&f1);
  
  QpGenerator qpGen;
  qpGen.init(&param, &oracle);
  
  while (!qpGen.run()) {

  }
  qpGen.output(cout);
  qpGen.terminate();
}

