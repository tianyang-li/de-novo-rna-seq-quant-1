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

using namespace Accpm;
using namespace std;

class LAOracleFunction : public OracleFunction {

private:
  AccpmGenMatrix *_A;
  AccpmVector *_b;

public:
  LAOracleFunction() : OracleFunction() {}
  LAOracleFunction(AccpmGenMatrix *A,  AccpmVector *b) : OracleFunction(), _A(A), _b(b) {}
  
  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {
   
    AccpmVector Ay(_b->size());
    AccpmLAMatVecMult(*_A, y, Ay);
    functionValue = 0.5 * AccpmLADotProd(Ay, y) - AccpmLADotProd(*_b, y);
   
    //  A * y - b;
    AccpmLAAddMult(Ay, -1, *_b); 
    memcpy(subGradients.addr(), Ay.addr(), sizeof(double)*Ay.size());

    if (info) {
      *info = 1;
    }
    return 0;
  }
};

static void createSDM(AccpmGenMatrix &A, AccpmVector &b);

int 
main(int argc, char *argv[])
{
  char *paramFile = "param.txt";
  bool addEqC = false;
  if (argc > 1) {
    if (argc != 2 && argc != 3) {
      std::cout << "Error starting " << argv[0] << std::endl;
      std::cout << "Usage: " << argv[0] << " <Parameter File> adeq[Add Equality Constraint]" 
	<< std::endl;
      exit(0);
    }
    paramFile = argv[1];
    if (argc == 3) {
      strcmp(argv[2], "adeq") == 0 ? addEqC = true : addEqC = false; 
    }
  } 
  
  Accpm::Parameters param(paramFile);
  int n = param.getIntParameter("NumVariables");
  
  vector<double> start(n, 0);
  /*
    if (n == 5) {
    start[0] = 5.0/6;
    start[1] = 2.0/3;
    start[2] = 0.5;
    start[3] = 1.0/3;
    start[4] = 1.0/6;	
    }
  */

  param.setStartingPoint(start);
  vector<double> center(n, 0.5);
  param.setCenterBall(center);

  if (addEqC) {
    AccpmGenMatrix constraint(n,2);
    constraint = 0;
    constraint(1,0) = 1;
    constraint(3,0) = 1;
    constraint(2,1) = 1;
    AccpmVector rhs(2);
    rhs(0) = 1;
    rhs(1) = 0.5;
    param.addEqualityConstraints(constraint, rhs);
    AccpmGenMatrix constraint2(n,1);
    constraint2 = 0;
    constraint2(0,0) = 1;
    constraint2(4,0) = 1;
    AccpmVector rhs2(1);
    rhs2 = 1;
    param.addEqualityConstraints(constraint2, rhs2);
  }

  AccpmGenMatrix A(n, n);
  AccpmVector b(n);
  // Initialize A as the symmetric difference matrix
  // and b as an identity vector
  createSDM(A, b);
  LAOracleFunction f1(&A, &b);
  Oracle oracle(&f1);

  QpGenerator qpGen;
  qpGen.init(&param, &oracle);
  while (!qpGen.run()) {
  }
  qpGen.output(cout);
  qpGen.terminate();
}

static
void createSDM(AccpmGenMatrix &A, AccpmVector &b)
{
  int n = b.size();
  A = 0;
  b = 0;
  b(0) = 1;
  A(0,0) = 2; 
  A(0,1) = -1;
  A(n-1,n-1) = 2;
  A(n-1,n-2) = -1;
  for(int i = 1; i < n-1; ++i) {

    A(i, i-1) = -1;
    A(i,i) = 2;
    A(i, i+1) = -1;
  }
  //cout << "A:\n" << A << endl;
  //cout << "b:\n" << b << endl;

}
