// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "AccpmVector.h"
#include "AccpmGenMatrix.h"
#include "AccpmDynMatrix.h"

#include <iostream>

using std::cout;
using std::endl;
static void test0LAVD();
static void test1LAGMD();
static void test2LAGMD(); // Test pre-allocation
static void test3LAGMD(); // Test memcpy vs inject
static void test4TRANSPOSE(); // Test memcpy vs inject
static void test5DYNMAT(); // Test dynamic matrix class

static int m = 10;
static int n = 3;
static bool verbose = true;
/*
static int m = 200;;
static int n = 10000;
static bool verbose = false;
*/

int 
main(int argc, char *argv[])
{
  int testNum = -1;
  if (argc > 1) {
    testNum = atoi(argv[1]);
  }
  switch (testNum) {
  case 0:
   test0LAVD();
   break;
  case 1:
    test1LAGMD();
    break;
  case 2:
    test2LAGMD();
    break;
  case 3:
    test3LAGMD();
    break;
  case 4:
    test4TRANSPOSE();
    break; 
  case 5:
    test5DYNMAT();
    break;
  default:
    test0LAVD();
    test4TRANSPOSE();
    test5DYNMAT();
  }
}

void 
test0LAVD()
{
  cout << "Test0: LAVD Vector Operations\n" << endl;
  Accpm::AccpmVector v1(1);
  v1 = 1;
  cout << "v1:\n" << v1 << endl;
 
  Accpm::AccpmVector v2(2);
  v2 = 2;
  cout << "v2:\n" << v2 << endl;
  v2.append(v1);
  cout << "v2.append(v1):\n" << v2 << endl;
  v1.append(v2);
  cout << "v1.append(v2):\n" << v1 << endl;  
  double v3 = 3;
  v1.append(v3);
  cout << "v1.append(v3):\n" << v1 << endl;  
  v1.negate();
  cout << "v1.negate():\n" << v1 << endl;
  v2.invert();
  cout << "v2.invert():\n" << v2 << endl;
  v2.times(v2);
  cout << "v2.times(v2):\n" << v2 << endl;
  v2.rdivide(v2);
  cout << "v2.rdivide(v2):\n" << v2 << endl;

  double sum = v1.sum();
  cout << "v1.sum():\n" << sum << endl;
  assert(sum == -9);
  double min = v1.min();
  cout << "v1.min():\n" << min << endl;
  assert(min == -3);
  double max = v1.max();
  cout << "v1.max():\n" << max << endl;
  assert(max == -1);
  if (v1 == v2) {
    cout << "v1 == v2" << endl;
    assert(0);
  } else {
    cout << "v1 != v2" << endl;
  }
  Accpm::AccpmVector v4(v1);
  if (v1 == v4) {
    cout << "v1 == v4" << endl;
  } else {
    cout << "v1 != v4" << endl;
    assert(0);
  }
  cout << "v4:\n" << v4 << endl;
  v4.deleteElem(0);
  v4.deleteElem(1);
  v4.deleteElem(v4.size()-1);
  cout << "Deleted v4(0), v4(1), v4(end):\n" << v4 << endl;
}

void
test1LAGMD()
{
  cout << "Test1: GMD Matrix\n" << endl;
  Accpm::AccpmVector ve(n);
  ve = 0;
  Accpm::AccpmGenMatrix A(ve);
  for (int i = 1; i < m; ++i) {
    Accpm::AccpmGenMatrix B(A.size(0),A.size(1)+1);
    B(A.index(0),A.index(1)).inject(A);
    ve = i;
    B.assignColumn(i, ve);
    if (verbose) {
      // cout << B.info() << endl;
      cout << "B:\n" << B << endl;
    }
    A = B;
  }
  cout << A.info() << endl << endl;
  for (int i = 0; i < A.size(1); ++i) {
    cout << "Column " << i << ":\n" << A.getColumn(i) << endl;
  }
  for (int i = 0; i < A.size(0); ++i) {
    cout << "Row " << i << ":\n" << A.getRow(i) << endl;
  }
  A.scaleColumn(9, 2);
  cout << "A: Column 9 scaled by 2\n" << A << endl;
  Accpm::AccpmGenMatrix B = A;
  A.addMult(2, B);
  cout << "A.addMult(2, A)\n" << A << endl;
  ve(0) = 2;
  ve(1) = 4;
  ve(2) = 6;
  cout << "ve.size()" << ve.size() << endl;
  A.scaleColumn(2, ve);
  cout << "A.scaleColumn(2, ve)\n" << A << endl;
  A.scale(3);
  cout << "A.scale(3)\n" << A << endl;
}


void
test2LAGMD()
{
  cout << "Test2: Resizing Memory pool\n" << endl;
  Accpm::AccpmVector ve(n);
  ve = 0;
  Accpm::AccpmGenMatrix A(ve);
  Accpm::AccpmGenMatrix B(A.size(0), m/2+1);
  B = -1;
  B(A.index(0),A.index(1)).inject(A);
  
  for (int i = 1; i < m; ++i) {
    
    if (i >= B.size(1)) {
      cout << "Reallocating" << endl;
      Accpm::AccpmGenMatrix newB(B.size(0),B.size(1)+m/2);
      newB = -1;
      newB(B.index(0),B.index(1)).inject(B);
      B = newB;
    }
    ve = i;
    B.assignColumn(i, ve);
    LaIndex colIndex(0, i);
    Accpm::AccpmGenMatrix subB = B(B.index(0), colIndex);
    if (verbose) {
      //cout << subB.info() << endl;
      cout << "subB:\n" << subB << endl;
      //cout << B.info() << endl;
      //cout << "B:\n" << B << endl;
    }
  }
  cout << B.info() << endl;
}

void
test3LAGMD()
{
  cout << "Test3: Memcpy\n" << endl;
  Accpm::AccpmVector ve(n);
  ve = 0;
  Accpm::AccpmGenMatrix *A = new Accpm::AccpmGenMatrix(ve);
  Accpm::AccpmGenMatrix *B;
  for (int i = 1; i < m; ++i) {
    B = new Accpm::AccpmGenMatrix(A->size(0),A->size(1)+1);
    //B(A.index(0),A.index(1)).inject(A);
    memcpy(B->addr(),A->addr(),sizeof(double)*A->size(0)*A->size(1));
    ve = i;
    B->assignColumn(i, ve);
    if (verbose) {
      // cout << B->info() << endl;
      cout << "B:\n" << *B << endl;
    }
    
    delete A;
    A = B; 
    /*
      int numRow = B.size(0);
      int numCol = B.size(1);
      A = new AccpmGenMatrix(numRow, numCol);
      memcpy(A->addr(),B.addr(),sizeof(double)*B.size(0)*B.size(1));
    */
  }
  cout << A->info() << endl << endl;
  delete A;
}

 
void test4TRANSPOSE()
{ 
  cout << "Test4: Transpose\n" << endl;
  Accpm::AccpmGenMatrix A(n, m);
  Accpm::AccpmVector ve(n);
  for (int i = 0; i < m; ++i) {
    ve = i;
    A.assignColumn(i, ve);
  }
  if (verbose) {
    cout << "A:\n" << A << endl;
  }
  cout << A.info() << endl << endl;
  Accpm::AccpmGenMatrix *AT = A.transpose();
  if (verbose) {
    cout << "AT:\n" << *AT << endl;
  }
  cout << AT->info() << endl << endl;
  Accpm::AccpmGenMatrix *ATT = AT->transpose();
  if (verbose) {
    cout << "ATT:\n" << *ATT << endl;
  }
  cout << ATT->info() << endl << endl;
  delete AT;
  delete ATT;
}

void
test5DYNMAT()
{
  cout << "Test5: Dynamic Matrix\n" << endl;
  Accpm::AccpmVector ve(n);
  ve = 0;
  Accpm::AccpmGenMatrix A(ve);
  Accpm::AccpmDynMatrix B1(false, true, A.size(0), 0);
  
  for (int i = 0; i < m; ++i) {
    ve = i;
    B1.addColumn(ve);
    Accpm::AccpmGenMatrix subB1 = B1.getM(); 

    if (verbose) {
      cout << subB1.info() << endl;
      cout << "subB1:\n" << subB1 << endl;
    }
  }
  B1.Info(cout);
  if (verbose) {
    cout << "B1: (Preallocated Columns)\n" << B1 << endl;
  }

  Accpm::AccpmDynMatrix B2 = B1;
  B2.Info(cout);
  if (verbose) {
    cout << "B2: (copy of B1)\n" << B2 << endl;
  }
  B2 = 1.7;
  if (verbose) {
    cout << "B1:\n" << B1 << endl;
    cout << "B2:\n" << B2 << endl;
  }

  Accpm::AccpmDynMatrix B3(true, false);
  for (int i = 0; i < m; ++i) {
    ve = i;
    B3.addRow(ve);
    Accpm::AccpmGenMatrix subB3 = B3.getM(); 

    if (verbose) {
      cout << subB3.info() << endl;
      cout << "subB3:\n" << subB3 << endl;
    }
  }
  B3.Info(cout);
  if (verbose) {
    cout << "B3: (Preallocated Rows)\n" << B3 << endl;
  } 
  B1.deleteColumn(2);
  cout << "B1.deleteColumn(2)\n" << B1 << endl;
  cout << "Deleted Column 2: new subB1:\n" << B1.getM() << endl;
  B1.deleteColumn(0);
  cout << "Deleted Column 1: new subB1:\n" << B1.getM() << endl;
  ve = 2;
  B1.addColumn(ve);
  cout << "Added Column: new subB1:\n" << B1.getM() << endl;
}
