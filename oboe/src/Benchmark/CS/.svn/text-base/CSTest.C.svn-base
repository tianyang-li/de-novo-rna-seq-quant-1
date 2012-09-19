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

static void
printUsage(const char *progName)
{
   std::cout << "Usage: " << progName << " <Parameter File> <Input Data File>"
	     << " <Number of variables>" 
	     << std::endl;
}

class CSOracleFunction : public OracleFunction {

private:
  AccpmVector _l;
  AccpmVector _lb;
  double _L;
 
  int readVector(std::ifstream &fin, AccpmVector &vector) {
    string line;
    std::istringstream instream;
    double val;
    int status = 0;
    for (int i = 0; i < vector.size(); ++i) {
      if (getline(fin, line)) {
	instream.clear();
	instream.str(line); 
	if ((instream >> val >> std::ws)) {
	  vector(i) = val; 
	} else {
	  status = 1;
	}
      } else {
	status = 1;
      }
    }
    return status;
  }

  void loadData(const string &fileName) {
    
    std::ifstream fin(fileName.c_str());
    if (fin.is_open()) {
      int status = 1;
      string line;
      std::istringstream instream;
      double val;
      int n;
      if (getline(fin, line)) {
	instream.clear();
	instream.str(line); 
	if ((instream >> n >> std::ws) && getline(fin, line)) {
	  assert(n == _l.size());
	  instream.clear();
	  instream.str(line); 
	  if ((instream >> val >> std::ws)) {
	    _L = val;
	    status = 0;
	  } 
	}
      }
      if (!status) {
	status |= readVector(fin, _l);
	status |= readVector(fin, _lb);
	_lb.negate();
      } 
      if (status) {
	std::cout << "\nError reading from file: " << fileName << std::endl;
	exit(1);
      }
    } else {
      std::cout << "\nError opening file: " << fileName << std::endl;
      printUsage("oboeCS");
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

  /**
   * Class to facilitate getting the indices of a vector of IndexedELement sorted on value. 
   */
  class IndexedElement {
  public:
    unsigned int index;
    double value;
  };

  static int
  cmpFunc(const void *a, const void *b)
  {
    IndexedElement *ia = (IndexedElement *)a;
    IndexedElement *ib = (IndexedElement *)b;
    if (ia->value == ib->value) {
      return 0;
    } 
    if (ia->value < ib->value) {
      return 1;
    } 
    return -1;
  }
  
  double
  findShortCut()
  {
    double amax = ACCPM_MINUS_INF;
    int index = -1;
    for (int i = 0; i < _l.size(); ++i) {
      if (amax < _l(i)) {
	amax = _l(i);
	index = i-1;
      } else {
	index = i;
      }
    }
    
    int index1 = -1;
    for (int i = 0; i < _l.size(); ++i) {
      if (amax > _l(i)) {
	index1 = i;
      }
    }
    assert(index >= 0);
    assert(index1 == index);
    double shortcut = (amax - 1)*(index + 1);
    return shortcut;

  }
  /**
   * KnapSacPart: 
   * Solves max c^Tz s.t. a^Tz <= lb
   **/
  void
  knapSacPart(const AccpmVector &a, const AccpmVector &c, int la, int lb, double shortcut,
	      AccpmVector &iz, AccpmVector &z)
  {
    double temp;
    double total;
    int t;
    for(int j = la-1; j < lb; ++j) {
      double maxvalue = 0;
      for(int i = 0; i < a.size(); ++i) {
	if((a(i)-1 <= j) && ( (i<=0) || (j-(a(i)-1) < shortcut))){
	  if(j-(a(i)-1) == 0) {
	    temp=0;
	  } else {
	    t = (int)a(i)-1;
	    temp = z(j-t -1);
	  }
	  total = c(i) + temp;
	  if (total >= maxvalue){
	    maxvalue = total;
	    iz(j) = i + 1;
	  }
	}
      }
      z(j) = maxvalue;
    }
  }

  void
  knapSac(const AccpmVector &y, AccpmVector &orig_x, double &obj)
  {
    int lb = int(_L);
    AccpmVector z(lb);
    AccpmVector iz(lb);
    z = 0;
    iz = 0;
    AccpmVector x(y.size());
    orig_x = AccpmVector(y.size());
    x = 0;
    IndexedElement ybyl[y.size()]; 
    for (int i = 0; i < y.size(); ++i) {
      ybyl[i].index  = i;
      ybyl[i].value  = y(i)/_l(i);
    }
   
    qsort(&ybyl, y.size(), sizeof(IndexedElement), cmpFunc);   
   
    AccpmVector c(y.size());
    AccpmVector a(y.size());
    for (int i = 0; i < y.size(); ++i) {
      c(i) = y(ybyl[i].index);
      a(i) = _l(ybyl[i].index);
    }
    double la;
    int index;
    findMin(a, la, index);
    double shortcut = findShortCut();
    knapSacPart(a, c, (int)la, lb, shortcut, iz, z);
    
    obj = z(lb-1);
    int k = (int)_L-1;
    double low = la - 1;
    while ((k+1 > low) && (iz(k) > 0)) {
      ++x((int)iz(k)-1);
      k = k - int(a((int)iz(k)-1));
    }
    for (int i = 0; i < y.size(); ++i) {
      orig_x(ybyl[i].index) = x(i);
    }
  }


public:
  CSOracleFunction() : OracleFunction() {}
  CSOracleFunction(const char *dataFile, int n) : OracleFunction(), _l(n), _lb(n) {
    loadData(dataFile);
  }

  
  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {
    
    double obj;
    AccpmVector x;
    knapSac(y, x, obj);
    if (obj <= 1) { /* Add optimality cut, no attractive pattern found for current dual variable y*/
      functionValue = AccpmLADotProd(_lb, y);
      memcpy(subGradients.addr(), _lb.addr(), sizeof(double)*_lb.size());
      *info = 1;
    } else { /* Add feasibility cut, which is the new pattern generated */
      functionValue = obj - 1;
      memcpy(subGradients.addr(), x.addr(), sizeof(double)*x.size());
      *info = 0;
    }

    return 0;
  }
};



int 
main(int argc, char *argv[])
{
  int n = 20;
  
  char *paramFile = "param.txt";
  char *dataFile = "CSData.dat";
  
  if (argc > 1) {
    if (argc != 4) {
      std::cout << "\nError starting " << argv[0] << std::endl;
      printUsage(argv[0]);
      exit(0);
    }
    paramFile = argv[1];
    dataFile = argv[2];
    n = atoi(argv[3]);
  } 
  Accpm::Parameters param(paramFile);
  if (n != param.getIntParameter("NumVariables")) {
    n = param.getIntParameter("NumVariables");
  }
  vector<double> start(n, 0.5);
  param.setStartingPoint(start);

  CSOracleFunction f1(dataFile, n);
  Oracle oracle(&f1);
      
  QpGenerator qpGen;
  qpGen.init(&param, &oracle);
      
  while (!qpGen.run()) {

  }
  qpGen.output(cout);
  cout << "The patterns are:" << endl;
  const AccpmGenMatrix *cuts;
  const AccpmVector *x;
  qpGen.getActiveCuts(cuts);
  x = qpGen.getCurrentX();
  cout << *cuts << endl;
  cout << "With weights:" << endl;
  cout << *x << endl;

  qpGen.terminate();
}

