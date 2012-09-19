// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef ORACLE_H
#define ORACLE_H

/**
 * Abstract class used for communicating to the query point
 " generator Accpm::QpGenerator. 
 *
 **/  

#include "AccpmDefs.h"
#include "AccpmGenMatrix.h"
/** 
 * @defgroup UserInterface
 * Group defining the essential classes which are used by the user 
 * directly. 
 * The user needs to use the Oracle, Parameters and QpGenerator objects in
 * order to implement the "Oracle" for the Oracle Based Optimization Engine
 * OBOE.
 *  
 */

/**
 * @file Oracle.h
 * @ingroup UserInterface
 */

namespace Accpm {

class OracleFunction {
  
 public:
  virtual ~OracleFunction(){};
  /**
   * The evaluate function which is called by OBOE.
   * At every iteration, i.e. call to QpGenerator::run() function,
   * OBOE asks the Oracle to provide it with the cut information
   * at the current query point y.
   * 
   * For a given query point, y, the eval function should return the following:
   *
   * @return functionValue which is an vector of size NumCuts x 1
   *    This value represents the following:
   *
   *    For the non-smooth function f1(.) :
   *    For optimality cuts, i.e *info = 1, functionValue is typically the
   *    function evaluation at the current query point, y.
   *    For feasiblity cuts, i.e *info = 0, however, this is a value which 
   *    would make the cut described by the subGradients matrix valid in the 
   *    form :
   *    functionValue + subGradients^T*(y' - y) <= 0, for all feasible points y'
   *
   *    For the smooth function f2(.) :
   *    It again represents the function value f2(y). 
   *
   * @return Sub-gradient at y in matrix subGradients which of size
   *    NumVariables x NumCuts. 
   *    The user can provide more than one cut for each query point, which
   *    usually is required when the NumSubProblems is more than 1.
   *
   *    The subgradient vector(with abuse of name) is also used to
   *    provide a valid cuts incase the given point, y, is not feasible.
   *    Hence the eval function is responsible for providing both optimality
   *    and feasibility cuts.
   *
   * @return AccpmGenMatrix *info:
   *    For non-smooth function f1 it is a vector of dimension NumCuts x 1.
   *    It specifies the cut type for each cut:
   *    (i) For feasibility cut,   *info(i,0) = 0 
   *    (ii) For optimality cut i, *info(i,0) = i, the index of subproblem.
   *    
   *	For the smooth function f2, this vector has the Hessian information. 
   *    It is a vector of dimension NumVariables x NumVariables 
   *    (if parameter diagHessian is false) otherwise its
   *	has dimension NumVariables x 1
   * 
   *	The return value is 0 on success and 1 to terminate the Outer Iterations
   *    and hence the Query point generation process.
   *
   *    Note: The number of cuts, NumCuts, returned at every iteration is 
   *    controlled by the user by providing the vectors and matrices of 
   * 	corresponding size.
   *    
  */
  virtual int eval(const AccpmVector &y, 
		   AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, 
		   AccpmGenMatrix *info) = 0;
};

/**
	* 	Oracle class is the container for the non-smooth, f1(.) and smooth function 
	*	f2(.), both of which are implemented via the interface OracleFunction.
*/

class Oracle {
 protected:
  OracleFunction *_f1;
  OracleFunction *_f2;

 public:
  /**
   * Constructor for Oracle takes in the non-smooth OracleFunction f1
   * and the smooth OracleFunction f2 if it is present.
   * The user needs to create their implementation of OracleFunction and
   * pass it to Oracle object.
   * If the application does not demand overriding the other virtual functions
   * of class Oracle, the user can use this class.
   * For cases when the user specifies a lower bound 
   * (for a minimization problem) they need to create their implementation of
   * Oracle class (in addition to OracleFunction which always needs to be
   * implemented by the user) and implement the computesLowerBound() 
   * and double getLowerBound() functions.
   *
   * @param f1 pointer to the implementation of user-defined OracleFunction
   *                      for the non-smooth part of the objective function  
   *
   * @param f2 pointer to the smooth part of the objective function.
  **/
  Oracle(OracleFunction *f1, OracleFunction *f2 = NULL);
  virtual ~Oracle();
  virtual bool hasSmoothOracle() const { return _f2 != NULL; }
  /**
   * Re-implement this function to return true give a user-defined lower bound
   * (upper bound) for minimization(maximization) problems.
   *
   **/
  virtual bool computesBound() const { return false; }
  /**
   * Re-implement this function to give a user-defined lower bound
   * upper bound for minimization(maximization) problems.
   *
   * Note: The function values returned in the OracleFunction::eval() are used to
   * determine the upper bound(lower bound) for minimization(maximization) problems.
  **/
  virtual double getObjectiveFunctionBound() const { return ACCPM_MINUS_INF; }
  
  /**
   * @return The non-smooth OracleFunction f1.
   **/
  OracleFunction *getF1() const { return _f1; }
   
  /**
   * @return The smooth OracleFunction f2.
   **/
  OracleFunction *getF2() const { return _f2; }

};
}
#endif
