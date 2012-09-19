// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <functional>
using std::string;
using std::vector;
using std::map;

#include "AccpmDefs.h"
#include "AccpmVector.h"
#include "AccpmGenMatrix.h"

/**
 *  Group for managing the user defined parameters used for 
 *  configuring OBOE.
 */
/**
 * @file Parameters.h
 * @ingroup UserInterface
 */

namespace Accpm 
{
  class Oracle;
  enum OptType { OPT_MIN = 1, OPT_MAX = -1 };
/**
 * Manages the parameters used for tuning OBOE. It can read
 * the parameters from a file with a specified format or can
 * be used to set the individual parameters.
 *
 * @ingroup ProblemInput */  

  class Parameters {
  
  private:   
    typedef bool (Parameters::*SetIntFunctor)(int);
    typedef map<const char *, SetIntFunctor, ltstr> SetIntFunctorMap;
    typedef int (Parameters::*GetIntFunctor)(void) const;
    typedef map<const char *, GetIntFunctor, ltstr> GetIntFunctorMap;
    
    typedef bool (Parameters::*SetRealFunctor)(Real);
    typedef map<const char *, SetRealFunctor, ltstr> SetRealFunctorMap;
    typedef Real (Parameters::*GetRealFunctor)(void) const;
    typedef map<const char *, GetRealFunctor, ltstr> GetRealFunctorMap;
 
    typedef bool (Parameters::*SetStringFunctor)(const string &);
    typedef map<const char *, SetStringFunctor, ltstr> SetStringFunctorMap;
    typedef const string &(Parameters::*GetStringFunctor)(void) const;
    typedef map<const char *, GetStringFunctor, ltstr> GetStringFunctorMap;

   
    string _problemName; 
    string _problemAcronym;
    string _optimizationType;
    OptType _type;
    unsigned int _numVariables;
    unsigned int _numSubProblems;
    
    Real _objLB;
    Real _objUB;

    AccpmVector *_startingPoint;
    AccpmVector *_varLB;
    AccpmVector *_varUB;
    AccpmVector *_b;
    AccpmVector *_pi;
    Real _delta;
    Real _eta;
    Real _gamma;
    Real _epsilonReal;
    Real _epsilonTol;

    string _methodName;
    int _maxOuterIterations;
    int _maxInnerIterations;
    int _verbosity;

    Real _tolerance;
    bool _filter;

    bool _convexityCheck;
    bool _convexityFix;
    bool _fixedProximalCenter;

    bool _proximal;
    Real _rho;
    bool _dynamicRho;
    Real _rhoMin;
    Real _rhoMax;
    Real _weightEpigraphCutInit;
    Real _weightEpigraphCutInc;
    
    bool _ball;
    Real _radiusBall;
    AccpmVector *_centerBall;
    
    Oracle *_oracle;
    bool _diagHessian;

    bool _computeLB;
    
    /** Verify the Localization Set has non-empty interior.
     * If the parameter is true, and localization set turns out to have
     * an empty interior, the program exits.
     */
    bool _checkLocSet;
    /**
     * The name of the lpsolver to use for checking the localization
     * set is non-empty. We use an OSI interface to the LP Solvers.
     **/
    string _lpSolverName;
    
    StdIntSet *_avIndex; // Active Variable Index
    
    AccpmGenMatrix *_D; // Equality Constraint
    AccpmVector *_d;    // rhs of Equality Constrain

    void init();
    void setFunctionMaps();

    SetIntFunctorMap _sintFunctionMap;
    GetIntFunctorMap _gintFunctionMap;
    SetRealFunctorMap _srealFunctionMap;
    GetRealFunctorMap _grealFunctionMap;
    SetStringFunctorMap _sstringFunctionMap;
    GetStringFunctorMap _gstringFunctionMap;

    // The int value access functions
    inline bool setNumVariables(int n) { _numVariables = n; return true; }
    inline int getNumVariables(void) const { return _numVariables; }
    inline bool setNumSubProblems(int n) { _numSubProblems = n; setPi(StdRealVector(n, 1)); return true; }
    inline int getNumSubProblems(void) const { return _numSubProblems; }
    inline bool setMaxOuterIterations(int n) { _maxOuterIterations = n; return true; }
    inline int getMaxOuterIterations() const { return _maxOuterIterations; }
    inline bool setMaxInnerIterations(int n) { _maxInnerIterations = n; return true; }
    inline int getMaxInnerIterations() const { return _maxInnerIterations; }
    inline bool setVerbosity(int n) { _verbosity = n; return true; }
    inline int getVerbosity() const { return _verbosity; }
    
    // Bool functions
    inline bool setFilter(int n) { n == 1 ? _filter = true : _filter = false;
      return true; }
    inline int getFilter() const { return _filter; }
    
    inline bool setConvexityCheck(int n) { n == 1 ? _convexityCheck = true 
	: _convexityCheck = false; return true; }
    inline int getConvexityCheck() const { return _convexityCheck; }
    inline bool setConvexityFix(int n) { n == 1 ? _convexityFix = true 
	: _convexityFix = false; return true; }
    inline int getConvexityFix() const { return _convexityFix; }
    inline bool setFixedProximalCenter(int n) { n == 1 ? _fixedProximalCenter = true 
	: _fixedProximalCenter = false; return true; }
    inline int getFixedProximalCenter() const { return _fixedProximalCenter; }

    bool setProximal(int n);
    int getProximal() const;

    inline bool setDynamicRho(int n) { n == 1 ? _dynamicRho = true 
	: _dynamicRho = false; return true; }
    inline int getDynamicRho() const { return _dynamicRho; }

    bool setBall(int n);
    inline int getBall() const { return _ball; }
    
    bool setBox(int n);
    inline int getBox() const { return _varLB || _varUB; }
    
    bool setDiagHessian(int n);
    inline int getDiagHessian() const { return _diagHessian; }
    
    inline bool setComputeLB(int n) { n == 1 ? _computeLB = true 
    	: _computeLB = false; return true; }
    inline int getComputeLB() const { return _computeLB; }
    
    inline bool setCheckLocSetInterior(int n) {  n == 1 ? _checkLocSet = true 
    	: _checkLocSet = false; return true; }
    inline int getCheckLocSetInterior() const { return _checkLocSet;}

    // The real value access functions
    inline bool setObjLB(Real value) { _objLB = value; return true; }
    inline Real getObjLB(void) const { return _objLB; }
    inline bool setObjUB(Real value) { _objUB = value; return true;}
    inline Real getObjUB(void) const { return _objUB; }
    inline bool setDelta(Real value) { _delta = value; return true;}
    inline Real getDelta(void) const { return _delta; }
    inline bool setEta(Real value) { _eta = value; return true;}
    inline Real getEta(void) const { return _eta; }
    inline bool setGamma(Real value) { _gamma = value; return true;}
    inline Real getGamma(void) const { return _gamma; }
    inline bool setTolerance(Real value) { _tolerance = value;return true; }
    inline Real getTolerance(void) const { return _tolerance; }
    inline bool setRho(Real value) { _rho = value; return true; }
    inline Real getRho(void) const { return _rho; }
    inline bool setRhoMax(Real value) { _rhoMax = value; return true; }
    inline Real getRhoMax(void) const { return _rhoMax; }
    inline bool setRhoMin(Real value) { _rhoMin = value; return true; }
    inline Real getRhoMin(void) const { return _rhoMin; }
    inline bool setWeightEpigraphCutInit(Real value) { _weightEpigraphCutInit = value; return true; }
    inline Real getWeightEpigraphCutInit(void) const { return _weightEpigraphCutInit; }
    inline bool setRadiusBall(Real value) { _radiusBall = value; return true; }
    inline Real getRadiusBall(void) const { return _radiusBall; }
    inline bool setEpsilonReal(Real value) { _epsilonReal = value; return true;}
    inline Real getEpsilonReal(void) const { return _epsilonReal; }
    inline bool setEpsilonTol(Real value) { _epsilonTol = value; return true;}
    inline Real getEpsilonTol(void) const { return _epsilonTol; }
    inline bool setWeightEpigraphCutInc(Real n) { _weightEpigraphCutInc = n; return true; }
    inline Real getWeightEpigraphCutInc() const { return _weightEpigraphCutInc; }
    
    inline bool setProblemName(const string &value) { _problemName = value; return true; }
    inline const string &getProblemName(void) const { return _problemName; }
    inline bool setProblemAcronym(const string &value) { _problemAcronym = value; return true; }
    inline const string &getProblemAcronym(void) const { return _problemAcronym; }
    bool setMethodName(const string &value);
    inline const string &getMethodName(void) const { return _methodName; }

    bool setLPSolverName(const string &value);
    inline const string &getLPSolverName(void) const { return _lpSolverName; }
   
    bool addParameter(const string &paramName, const string &paramType, const string &paramValue);
   
  public:
  
    Parameters(const char *fileName);
    virtual ~Parameters();
    
    bool setIntParameter(const char *name, int value);
    int getIntParameter(const char *name) const;
    bool setRealParameter(const char *name, Real value);
    Real getRealParameter(const char *name) const;
    bool setStringParameter(const char *name, const string &value);
    const string getStringParameter(const char *name) const;

    bool setOptimizationType(const string &value);
    const string &getOptimizationType(void) const;
    const OptType getOptType(void) const { return _type; }

    /**
     * Set the starting point for the problem.
     * This sets the first query point given to the Oracle::eval() function
     * by OBOE.
     * It is the user responsibility to ensure the starting point is feasible
     * with respect to the variable bounds which are set by setVariableLB()
     * and setVariableUB().
     *
     * @return Returns true if it is successful in setting the specified point.
     *
     * If no starting point is specified the starting point is set to 0 for 
     * all the variables.
     **/
    bool setStartingPoint(const StdRealVector& v);

    /**
     * Set the lower bound on the variables.
     * Currently the user needs to specify the bound on all the variables if
     * they choose to set the lower bound. Ofcourse the bound can be set 
     * ACCPM_MINUS_INF for variables with no lower bounds.
     * 
     * The default lower bound on all variables is ACCPM_MINUS_INF.
     **/
    bool setVariableLB(const StdRealVector& v);

    /**
     * Set the upper bound on the variables.
     * Currently the user needs to specify the bound on all the variables if
     * they choose to set the upper bound. Ofcourse the bound can be set 
     * ACCPM_PLUS_INF for variables with no upper bounds.
     * 
     * The default lower bound on all variables is ACCPM_PLUS_INF.
     **/
    bool setVariableUB(const StdRealVector& v);
    
    /**
     * Set the vector of weights on the non-smooth function, f1.
     * The vector pi should be the same size as the NumSubProblems.
     **/
    bool setPi(const StdRealVector& pi);
    
    /**
     * Set the linear component of the objective function.
     * Note: The user can either specify the linear component or
     * provide a smooth oracle, f2, but not both.
     * If the smooth function f2 is completely linear then the user
     * should use the setB() function and not write a smooth oracle.
    **/
    bool setB(const StdRealVector& v);
    
    /**
     * Same as setB(const StdRealVector& v), with the difference that
     * the linear component is specified as an AccpmVector instead of
     * vector<double>.
     **/
    bool setB(const AccpmVector& v);
    
    /**
     * Set the center of the ball, yr, if there are ball constraints:
     *                   |y - yr|^2 <= R^2
     * where R is specified by RadiusBall parameter.
    **/
    bool setCenterBall(const StdRealVector& v);
  
    inline const AccpmVector *getStartingPoint() const { return _startingPoint; }
    inline const AccpmVector *getVariableLB() const { return _varLB; }
    inline const AccpmVector *getVariableUB() const { return _varUB; }
    inline const AccpmVector *getPi() const { return _pi; }
    inline const AccpmVector *getB() const { return _b; }
    inline const AccpmVector *getCenterBall() const { return _centerBall; }
    
    bool setOracle(Oracle *oracle);
    const Oracle *getOracle(void) const; 

    /**
     * Set the variables whose indices are specified in the set v to be active,
     * all the other variables are considered to be put to 0 and hence
     * de-activated.
     **/
    bool setActiveVariableIndex(const StdIntSet &v);
    /**
     * @return The index of the active variables at the current iteration.
     */
    const StdIntSet *getActiveVariableIndex() const;
    /**
     * @return True if the dimension of the problem can change from one iteration to the
     * next. This is possible if some of the variables are de-activated by not making them
     * part of the active set specified by Parameters::setActiveVariableIndex();
     **/
    bool hasVariableDimension() const;
 
    /**
     * Add Equality constraints (D^t)*y = d
     * The constraints matrix represents D and should have dimension 
     * number of variables(n) * number of constraints.
     * The vector rhs is the vector with dimension number of constraints and 
     * represents the vector d, right hand side of the constraint.
     **/
    void addEqualityConstraints(const AccpmGenMatrix &constraints, const AccpmVector &rhs);
    /**
     * Remove the constraints (D^t)*y = d if they exist
     */
    void removeEqualityConstraints();
    /**
     * Whether the problem has Equality Constraints of the type (D^t)*y = d
     */
    bool hasEqualityConstraints() const;
    /**
     * @return A pointer to Equality constraint matrix D^t and corresponding rhs vector d
     * if constraints (D^t)*y = d have been added.
     */
    void getEqualityConstraints(AccpmGenMatrix *&constraints, AccpmVector *&rhs) const;

    void output(std::ostream &os) const;
    //* I/O
	friend std::ostream& operator<<(std::ostream&, const Parameters &P);
  };
}

#endif
