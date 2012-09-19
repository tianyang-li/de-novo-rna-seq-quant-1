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
#include "QPNDHSmoothOF.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace Accpm;
using namespace std;
using std::string;

class QPOracleFunction : public OracleFunction {

public:
  AccpmGenMatrix _A;
  AccpmGenMatrix _B;
  AccpmGenMatrix _H;
  AccpmGenMatrix _F;
  int _formulation;

private:
  void readMatrix(const string &fileName, AccpmGenMatrix &matrix) {
    //std::cout << "Reading File:" << fileName << std::endl;
    std::ifstream fin(fileName.c_str());
    if (fin.is_open()) {
      string line;
      std::istringstream instream;
      double val;
      int n;
      int m;
      // Read the dimesion of the matrix;
      if (getline(fin, line)) {
	  instream.clear();
	  instream.str(line); 
	  instream >> n >> std::ws >> m;
      }
      matrix.resize(n, m);
      std::cout << "Reading matrix of dimesion: " << n << " * " << m << std::endl; 
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

  int eval1(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {
    int n = y.size();
    bool feasible = true;
   
    AccpmVector Hy(n);
    AccpmLAMatVecMult(_H, y, Hy);
    std::cout << "Hy: " << Hy << std::endl;
    std::cout << "F: " << _F << std::endl;
    // 1/2(y^THy) - d 
    double constraintGap = 0.5*AccpmLADotProd(y, Hy);
    constraintGap -= _F(0,0);
    std::cout << "gap: " << constraintGap << std::endl;
    if (constraintGap > 0) {
      feasible = false;
      
      memcpy(subGradients.addr(), Hy.addr(), sizeof(double)*n);
      functionValue = constraintGap;
      *info = 0;
    } else {
      AccpmVector Ay(_A.size(0));
      AccpmLAMatVecMult(_A, y, Ay);
      
      AccpmVector grad = _B.getColumn(0);
      functionValue = 0.5*AccpmLADotProd(Ay, Ay) + AccpmLADotProd(grad, y) ;

      AccpmLAMatTransVecMult(_A, Ay, grad, 1, 1);
      memcpy(subGradients.addr(), grad.addr(), sizeof(double)*n);

      *info = 1;
    }
    return 0;
  }

  int eval2(const AccpmVector &y, AccpmVector &functionValue, 
	    AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {
  
    int nx = _A.size(0);
    int ny = _A.size(1);
    AccpmVector yvar = y(LaIndex(nx, y.size()-1));
    AccpmVector Hy(ny);
    AccpmLAMatVecMult(_H, yvar, Hy);
    std::cout << "Hy: " << Hy << std::endl;
    std::cout << "F: " << _F << std::endl;
    // 1/2(y^THy) - d 
    double constraintGap = 0.5*AccpmLADotProd(yvar, Hy);
    constraintGap -= _F(0,0);
    if (DBL_GT(constraintGap, 0)) { //infeasible
      std::cout << "gap: " << constraintGap << std::endl;
      subGradients = 0;
      memcpy(subGradients.addr()+nx, Hy.addr(), sizeof(double)*ny);
      functionValue = constraintGap;
      *info = 0;
    } else {
      AccpmVector grad = _B.getColumn(0);
      functionValue = AccpmLADotProd(grad, yvar) ;
      subGradients = 0;
      memcpy(subGradients.addr()+nx, grad.addr(), sizeof(double)*ny);

      *info = 1;
    }
   
    return 0;
  }

public:
  QPOracleFunction() : OracleFunction() {}
  QPOracleFunction(const char *filePrefix, int formulation = 1) 
    : OracleFunction(), _formulation(formulation)
  {    
    string fileName(filePrefix);
    
    fileName = fileName + "A" +".dat";
    readMatrix(fileName, _A);
    
    fileName = string(filePrefix) + "B" + ".dat";
    readMatrix(fileName, _B);
    
    fileName = string(filePrefix) + "H" + ".dat";
    readMatrix(fileName, _H);
    
    fileName = string(filePrefix) + "F" + ".dat";
    readMatrix(fileName, _F);
  }
  
  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {
    if (_formulation == 1) {
      return eval1(y, functionValue, subGradients, info);
    } 
    if (_formulation == 2) {
      return eval2(y, functionValue, subGradients, info);
    } 
    std::cerr << "Unknown formulation" << std::endl;
    exit(0);
    return 0;
  } 
  
  
};

int 
main(int argc, char *argv[])
{
  char *paramFile = "paramNDH.txt";
  char *matrixFilePrefix = "QPNDH_";
  int formulation = 1;
  if (argc > 1) {
    if (argc != 3 && argc != 4) {
      std::cout << "Error starting " << argv[0] << std::endl;
      std::cout << "Usage: " << argv[0] << " <Parameter File> <Input Data File prefix>"
		<< " [Formulation][1,2] " << std::endl;
      exit(0);
    }
    paramFile = argv[1];
    matrixFilePrefix = argv[2];
    if (argc == 4) {
      formulation = atoi(argv[3]);
    }
  } 
  
  QPOracleFunction f1(matrixFilePrefix, formulation);
  int n;
  QPNDHSmoothOF *f2 = 0;
  int nx = f1._A.size(0);
  int ny = f1._A.size(1);
  if (formulation == 1) {
    n = ny;
  } else {
    n = nx + ny;
    f2 = new QPNDHSmoothOF(nx, ny);
  }
  Oracle oracle(&f1, f2);
  Accpm::Parameters param(paramFile);
  n = param.setIntParameter("NumVariables", n);

 /* 
  vector<double> start;
  start.push_back(-2.6059);
  start.push_back(-0.0081);
  start.push_back(-3.8686);
  param.setStartingPoint(start);
*/ 
  vector<double> varUB(n, 0);
  varUB[0] = 10;
  param.setVariableUB(varUB);
  vector<double> varLB(n, -10);
  param.setVariableLB(varLB);

  if (formulation == 2) {
    AccpmVector rhs(nx);
    rhs = 0;
    AccpmGenMatrix constraints(nx + ny, nx);
    constraints = 0;
    for (int i = 0; i < nx; ++i) {
      constraints(i, i) = -1;
    }
    for (int i = nx; i < nx + ny; ++i) {
      constraints.assignRow(i, f1._A.getColumn(i - nx));
    }
    std::cout << "Constraints: " << constraints << std::endl;
    param.addEqualityConstraints(constraints, rhs);
  }
  QpGenerator qpGen;
  qpGen.init(&param, &oracle);
  
  while (!qpGen.run()) {

  }
  qpGen.output(cout);
  std::cout << "Optimal:\n" << std::setprecision(16) << qpGen.getOptimalObj() 
	 << std::endl;
  std::cout << "Solution:\n" << *qpGen.getQueryPoint() << std::endl;
  qpGen.terminate();
  
  if (f2) {
    delete f2;
  }
}

