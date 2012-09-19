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

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace Accpm;
using namespace std;
using std::string;

class TSPOracleFunction : public OracleFunction {

private:
  AccpmGenMatrix _C;
  bool _coordInfo;

  void readMatrix(const string &fileName, AccpmGenMatrix &matrix) {
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
  
  double 
  rdist(int i, int j, const AccpmVector &y) {
    double rdist = 0;
    if (_coordInfo) {
      rdist = y(i) + y(j)
      + int(sqrt( pow(_C(i,0)-_C(j,0), 2)+ pow(_C(i,1)-_C(j,1),2))+0.5);
    } else {
      rdist = y(i) + y(j) + _C(i, j);
    }
    return rdist;
  }
  void
  findMin(const std::vector<double> &v, double &nmin, int &index) 
  {
    nmin = ACCPM_PLUS_INF;
    index = -1;
    for (unsigned int i = 0; i < v.size(); ++i) {
      if (nmin > v[i]) {
	nmin = v[i];
	index = i;
      }
    }
  }

  void 
  tspsim(const AccpmVector &y, double &fVal, AccpmVector &subGrad)
  {
    int n = y.size();
    std::vector<int> inxt(n, 0);
    std::vector<int> outxt(n, 0);
    std::vector<int> deg(n, 0);
    std::vector<double> tdist(n, 0);
    
    subGrad = AccpmVector(n);
    subGrad = 0;

    int k = 0;
    int outnode = n-1;
    for (int i = 0; i < n; ++i) {
      if (i != 0) {
	inxt[k] = 1;
	outxt[k] = i+1;
	tdist[k] = rdist(i, 0, y);
	++k;
      }
    }
    tdist.pop_back();
    double nmin;
    int index;
    findMin(tdist, nmin, index);
    assert(index >= 0);
    int ind = outxt[index];
    deg[ind-1] = 1;
    --outnode;
    double weight = nmin;
    int afterRoot = ind;
    
    k = 0;
    for (int i = 0; i < n; ++i) {
      if (i != 0 && i != ind-1) {
	inxt[k] = ind;
	outxt[k] = i+1;
	tdist[k] = rdist(i, ind-1, y);
	++k;
      }
    }

    while (outnode > 0) {
      tdist.pop_back();
      findMin(tdist, nmin, index);
      weight += nmin;
      int t1 = inxt[index];
      int t2 = outxt[index];
      inxt[index] = inxt[outnode-1];
      outxt[index] = outxt[outnode-1];
      tdist[index] = tdist[outnode-1];
      deg[t1-1] = deg[t1-1]+1;
      deg[t2-1] = deg[t2-1]+1;
      --outnode;
      for (int i = 0; i < outnode; ++i) {
	double e = rdist(t2-1,outxt[i]-1,y);
	if ( e < tdist[i]) {
	  inxt[i] = t2;
	  tdist[i] = e;
	}
      }
    }
    double min1 = ACCPM_PLUS_INF;
    int minin1 = -1;
    for (int i = 0; i < n; ++i) {
      if ( i!=0 && i != afterRoot-1) {
        double e = rdist(i, 0, y);
        if (e < min1) {
	  min1 = e;
	  minin1 = i;
	}
      }
    }
    weight += min1;
    fVal = -weight + 2*y.sum();
    
    deg[minin1] = deg[minin1]+1;
    deg[0] = 2;
    for (int i = 0; i < n; ++i) {
      subGrad(i) = 2 - deg[i];
    }
  }

public:
  TSPOracleFunction() : OracleFunction() {}
  TSPOracleFunction(const char *fileName, int n, bool coord) : OracleFunction(), 
							       _coordInfo(coord) {
    if (_coordInfo) {
      _C = AccpmGenMatrix(n, 2);
    } else {
      std::cout << "Using distance format" << endl;
      _C = AccpmGenMatrix(n, n);
    }
    readMatrix(fileName, _C);
  }

 

  virtual int 
  eval(const AccpmVector &y, AccpmVector &functionValue, 
       AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {
    
    double fVal;
    AccpmVector subGrad;

    tspsim(y, fVal, subGrad);

    functionValue = -fVal;

    subGrad.negate();
    memcpy(subGradients.addr(), subGrad.addr(), sizeof(double)*subGrad.size());

    *info = 1;

    return 0;
  }
};

int 
main(int argc, char *argv[])
{
  int n = 29;
  char *paramFile = "param.txt";
  string distanceFile;
  bool coord = true;

  if (argc > 1) {
    if (argc != 3 && argc != 4 && argc != 5) {
      std::cout << "Error starting " << argv[0] << std::endl;
      std::cout << "Usage: " << argv[0] << " <Parameter File> <Number of variables> <Distance File> [coord]"
		<< std::endl;
      exit(0);
    }
    paramFile = argv[1];
    n = atoi(argv[2]);
    if (argc > 3) {
      distanceFile = argv[3];
    }
    if (argc == 5) {
      coord = strcmp("coord",argv[4]) == 0;
    }
  } else {
    distanceFile = "tsp";
    ostringstream oss;
    oss << n;
    distanceFile += oss.str() + ".coord";
  }
  Accpm::Parameters param(paramFile);
  assert(n == param.getIntParameter("NumVariables"));
  vector<double> start(n, 0);
  param.setStartingPoint(start);
  param.setOptimizationType("Max");

  TSPOracleFunction f1(distanceFile.c_str(), n, coord);
  Oracle oracle(&f1);
 
    
 QpGenerator qpGen;
 qpGen.init(&param, &oracle);
  
 while (!qpGen.run()) {

  }
  qpGen.output(cout);
  qpGen.terminate();
}

