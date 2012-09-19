// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "blas1pp.h"
#include "blas2pp.h"
#include "Oracle.h"
#include "QpGenerator.h"
#include "Parameters.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include "EigenValue.h"
#include "laslv.h"

using namespace Accpm;
using namespace std;
using std::string;

/* Cant seem to resolve the conflict with the General Matrix equivalent.
 * Declaring here explicitly.
 **/
extern void Blas_Mat_Vec_Mult(const LaGenMatComplex &A, 
			      const LaVectorComplex &dx, 
			      LaVectorComplex &dy, 
			      LaComplex alpha = 1.0, LaComplex beta = 0.0);
extern void Blas_Mat_Trans_Vec_Mult(const LaGenMatComplex &A, 
				    const LaVectorComplex &dx, 
				    LaVectorComplex &dy,
				    LaComplex alpha = 1.0, LaComplex beta = 0.0);
extern double Blas_Norm2(const LaVectorComplex &dx);
extern void Blas_Scale(COMPLEX da, LaVectorComplex &dx);
extern COMPLEX Blas_H_Dot_Prod(const LaVectorComplex &cx, const LaVectorComplex &cy);
extern void Blas_R1_Update(LaGenMatComplex &A, const LaVectorComplex &dx, 
			   const LaVectorComplex &dy, LaComplex alpha = 1.0);

class NNPOLOracleFunction : public OracleFunction {

  typedef std::complex<double> dcomplex;

private:
  AccpmComplexMatrix _A;
  AccpmGenMatrix _rA;
  AccpmGenMatrix _iA;
  AccpmGenMatrix _B;
  AccpmGenMatrix _C;
  AccpmComplexVector _q;
  AccpmVector _minusB;
  
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

  void toeplitz(const AccpmComplexVector &s, AccpmComplexMatrix &M) 
  {
    for (int i = 0; i < s.size(); ++i) {
      for (int j = i; j < s.size(); ++j) {
	M(i,j) = s(j-i);
      }
    }
    for (int i = 0; i < s.size(); ++i) {
      for (int j = 0; j < i; ++j) {
	M(i,j) = LaComplex(LaComplex(M(j,i)).real(), - LaComplex(M(j,i)).imag());
      }
    }
  }
  /**
   * Dual Topelitz function
   */
  void dtoeplitz(AccpmComplexMatrix &M, AccpmComplexVector &f)
  {
    int n = f.size();
    f = 0;
    //for ind=1:n, I=ind+[0:n-ind]*(n+1); f(ind)=sum(Y(I)); end
    
    for (int i = 0; i < n; ++i) {
      AccpmVector I;
      for (int k = 0; k < n-i; ++k) {
	I.append(i + k*(n+1));
      }
      
      for (int j = 0; j < I.size(); ++j) {
	LaComplex c = LaComplex(*(M.addr() + int(I(j))));
	f(i) = LaComplex(f(i)) + c;
	// f(2:n)=2*f(2:n);
      }
      if (i != 0) {
	f(i) = LaComplex(2, 0) * LaComplex(f(i));
      }
    }
  } 

public:
  NNPOLOracleFunction() : OracleFunction() {}
  NNPOLOracleFunction(const char *filePrefix, int m, int n) : OracleFunction(), _A(m,n), _rA(m,n), _iA(m,n),
							      _B(m,1), _C(n,1), _q(n) {
    string fileName(filePrefix);
    fileName = string(filePrefix) + "rA" + ".dat";
    readMatrix(fileName, _rA);
    fileName = string(filePrefix) + "iA" + ".dat";
    readMatrix(fileName, _iA);
    
    fileName = string(filePrefix) + "B" + ".dat";
    readMatrix(fileName, _B);
    fileName = string(filePrefix) + "C" + ".dat";
    readMatrix(fileName, _C);

    /* Read the Q vector in order to have the same behaviour as versionM test */
    AccpmGenMatrix q(n, 1);
    fileName = string(filePrefix) + "Q" + ".dat";
    readMatrix(fileName, q);

    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
	_A(i,j) = LaComplex(_rA(i,j), _iA(i,j));
      }
    }
    for (int i = 0; i < _q.size(); ++i) {
      _q(i) = LaComplex(q(i, 0), 0);
    }
    /* for (int i = 0; i < n; ++i) {
      _q(i) = LaComplex(drand48(), 0);
      }*/
  
    _minusB = _B.getColumn(0);
    _minusB.negate();
    
  }
  
  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {
    
  
    //C - (real(A)' - i*imag(A)')*y
    LaVectorComplex s(_C.size(0));
    for (int i = 0; i < s.size(); ++i) {
      s(i) = LaComplex(_C(i,0), 0);
    }
    LaVectorComplex cy(y.size());
    for (int i = 0; i < y.size(); ++i) {
      cy(i) = LaComplex(y(i), 0);
    }

    Blas_Mat_Trans_Vec_Mult(_A, cy, s, -1, 1);
    AccpmComplexMatrix M(s.size(), s.size());
    toeplitz(s, M);
    
    AccpmVector eigenValues(s.size());
    if (computeEigenValues(M, eigenValues)) {
      std::cout << "Error computing Eigen Values of M." << std::endl;
    }
    double l = eigenValues.min();
    if (l >= 0) { // OptimalityCut
      functionValue = Blas_Dot_Prod(_minusB, y);
      subGradients = _minusB;
      
      *info = 1;
    } else {
      functionValue = -l;
      AccpmComplexVector q(_q);
      AccpmComplexMatrix Ms(M);
      for (int i = 0; i < M.size(0); Ms(i,i) = LaComplex(M(i,i)) - l, ++i);
      
      AccpmComplexVector z(Ms.size(1));
      while (1) {
	LaLinearSolve(Ms, z, q);
	q = z;
	Blas_Scale(LaComplex(1.0/Blas_Norm2(z), 0), q);
	AccpmComplexVector Mq(M.size(0)); 
	Blas_Mat_Vec_Mult(M, q, Mq, 1, 0);
	LaComplex lest = Blas_H_Dot_Prod(q, Mq);
	if (std::abs(l-lest.real()) < 1e-8) {
	  break;
	}
      }
      
      AccpmComplexMatrix qqT(q.size(), q.size());
      qqT = 0;
      Blas_R1_Update(qqT, q, q, 1);
      AccpmComplexVector f(q.size());
      dtoeplitz(qqT, f);
            
      AccpmVector subGrad(y.size());
      AccpmVector tmpf(f.size());
      for (int i = 0; i < tmpf.size(); ++i) {
	tmpf(i) = LaComplex(f(i)).real();
      }
      Blas_Mat_Vec_Mult(_rA, tmpf, subGrad, 1, 0);
      for (int i = 0; i < tmpf.size(); ++i) {
	tmpf(i) = LaComplex(f(i)).imag();
      }
      Blas_Mat_Vec_Mult(_iA, tmpf, subGrad, 1, 1);
      
      for (int i = 0 ; i < subGrad.size(); ++i) {
	subGradients(i, 0) = LaComplex(subGrad(i)).real();
      }
      *info = 0;
    }
    
    return 0;
  }
};

int 
main(int argc, char *argv[])
{
  int n = 26;
  int m = 10;
  char *paramFile = "param.txt";
  char *matrixFilePrefix = "NNPOL_";
  if (argc > 1) {
    if (argc != 5) {
      std::cout << "Error starting " << argv[0] << std::endl;
      std::cout << "Usage: " << argv[0] << " <Parameter File> <Input Data File prefix>"
		<< " <Number of variables> <Number of constraints>" 
		<< std::endl;
      std::cout << argv[0] << " needs to read the A, B, C matrices and Q vector" << std::endl;
      std::cout << "A is read in 2 matrices real part from <prefix>rA.dat and imaginary from"
		<< " <prefix>iA.dat" << std::endl;
      std::cout << "The other matrices are read as <prefix><B/C/Q>.dat" << std::endl;
      exit(0);
    }
    paramFile = argv[1];
    matrixFilePrefix = argv[2];
    m = atoi(argv[3]);
    n = atoi(argv[4]);
  } 
  Accpm::Parameters param(paramFile);
  assert(m == param.getIntParameter("NumVariables"));
  vector<double> start(m, 0);
  param.setStartingPoint(start);
  vector<double> varLB(m, -50);
  param.setVariableLB(varLB);
  vector<double> varUB(m, 50);
  param.setVariableUB(varUB);
  
  NNPOLOracleFunction f1(matrixFilePrefix, m, n);
  Oracle oracle(&f1);
  
  
 QpGenerator qpGen;
 qpGen.init(&param, &oracle);
 
  while (!qpGen.run()) {
    
  }
  qpGen.output(cout);
  qpGen.terminate();
  
  return 0;
}

