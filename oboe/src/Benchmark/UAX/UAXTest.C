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

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace Accpm;
using namespace std;
using std::string;

class UAXOracleFunction : public OracleFunction {

private:
  AccpmGenMatrix _A;
  AccpmGenMatrix _B;
 
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
  void
  findMin(const AccpmVector &v, double &nmin, int &index) 
  {
    nmin = ACCPM_PLUS_INF;
    index = -1;
    for (int i = 0; i < v.size(); ++i) {
      if (nmin > v(i)) {
	nmin = v(i);
	index = i;
      }
    }
  }

public:
  UAXOracleFunction() : OracleFunction() {}
  UAXOracleFunction(const char *filePrefix, int n, int m) : OracleFunction(), _A(n,1), _B(m-1,n) {
    
    string fileName(filePrefix);
    fileName = fileName + "A" +".dat";
    readMatrix(fileName, _A);
    fileName = string(filePrefix) + "B" + ".dat";
    readMatrix(fileName, _B);
  }
  
  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {
    AccpmVector c = _A;
    AccpmLAMatTransVecMult(_B, y, c, 1, 1);
    double cMin;
    int index;
    findMin(c, cMin, index);
    double gamma0 = y.sum() - 1;
    AccpmVector subGrad(y.size());

    if (gamma0 <= 0) {
      functionValue = c(index);
      subGrad = _B.getColumn(index);
      *info = 1; // optimality cut
    } else {
      functionValue = gamma0;
      subGrad = 1;
      *info = 0; // fesibility cut;
    }
    memcpy(subGradients.addr(), subGrad.addr(), sizeof(double)*subGrad.size());

    return 0;
  }
};

int 
main(int argc, char *argv[])
{
  int n = 10;
  int m = 10;
  char *paramFile = "param.txt";
  char *matrixFilePrefix = "UAX_";
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
  assert(m-1 == param.getIntParameter("NumVariables"));
  
  vector<double> start(m-1, 1.0/m);
  param.setStartingPoint(start);
  vector<double> varLB(m-1, 0);
  param.setVariableLB(varLB);
  vector<double> varUB(m-1, 1);
  param.setVariableUB(varUB);
  param.setOptimizationType("Max");
  
  UAXOracleFunction f1(matrixFilePrefix, n, m);
  Oracle oracle(&f1);
  
  QpGenerator qpGen;
  qpGen.init(&param, &oracle);
  
  while (!qpGen.run()) {

  }
  qpGen.output(cout);
  qpGen.terminate();
}

