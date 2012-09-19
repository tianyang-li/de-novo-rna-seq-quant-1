// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef METHOD_H
#define METHOD_H

/**
 * @file Method.h
 * @ingroup AccpmCore
 */
#include "AccpmVector.h"
#include "AccpmGenMatrix.h"
#include "Solution.h"
#include "symd.h"
#include "Parameters.h"

namespace Accpm 
{
  class Manager;
  class LocSet;
/**
 * Abstract class for methods generating query points.
 *
 * @ingroup AccpmCore */  

  class Method {
  
  protected:
    int _numIter;
    bool _computeLB;
    double _ratioCutVsDim;
    double _startLB;

  public:
    Method() : _numIter(0), _computeLB(true), _ratioCutVsDim(1), _startLB(1e-2) {};
    virtual int run(Manager &manager) = 0;
    virtual ~Method() {};
  };

  class DualMethod : public Method {
  
  private:

    const Parameters *_param;
    
    int performNewtonSteps(Manager &manager, LocSet &locSet);
    double computeLowerBound(Manager &manager, const LocSet &locSet);

    int computeInfeasibilities(const Manager &manager, const LocSet &locSet, 
			       const AccpmVector &y, const AccpmVector &z, double zs,
			       AccpmVector &s, double &s0,
			       AccpmVector &rs, double &rs0);
    void updateRsS(const Manager &manager, const LocSet &locSet,
		   const AccpmVector &ATy, const AccpmVector &z,
		   AccpmVector &rs, AccpmVector &s);
  
    int newtonDualDir(Manager &manager, LocSet &locSet, 
		      const AccpmVector &ry, const AccpmVector &rz,
		      const AccpmVector &rs, double rs0, 
		      const AccpmVector &wsV, double ws0, 
		      double f2, const AccpmVector *df2,
		      const AccpmGenMatrix *d2f2, 
		      const AccpmVector &diagQ,
		      const NewtonSolution &sol,
		      NewtonSolution &step, bool eqConstraint = false) const;
    
    int newtonDualDir1(Manager &manager, LocSet &locSet, 
		       const AccpmVector &ry, 
		       const AccpmVector &rs, double rs0,
		       const AccpmVector &rz, double rzs, double rss,
		       const AccpmVector &diagQ,
		       const AccpmVector *df2,
		       const NewtonSolution &sol,
		       NewtonSolution &step,
		       bool eqConstraint) const;

    int newtonDualDir2(Manager &manager, LocSet &locSet, 
		       const AccpmVector &ry, 
		       const AccpmVector &rs, double rs0,
		       const AccpmVector &wsV,
		       const AccpmVector &rz, double rzs, double rss,
		       const AccpmVector &diagQ,
		       const AccpmVector *df2,
		       const NewtonSolution &sol,
		       NewtonSolution &step) const;
    /**
     * Handled different factorization of the Newton System in case of large p.
     * Makes use of the Rank-1 updates.
     */
    int newtonDualDirForLargeP(Manager &manager, LocSet &locSet,
			       const SymmetricMatrix &H, const AccpmGenMatrix &H21, 
			       const AccpmVector &H13, const AccpmVector &ELET, double wsByssSq, 
			       double w0Bys0Sq, const AccpmVector &rhs1, const AccpmVector &rhs2, double rhs3, 
			       NewtonSolution &step) const;

    int computeRy(const Manager &manager, const LocSet &locSet, const AccpmVector &y,
		  const AccpmVector &wsV, const AccpmVector *df2, double ss, AccpmVector &ry) const;
    int computeRz(const LocSet &locSet, const AccpmVector &wsV, double ws0, AccpmVector &rz) const;
   
    int prodAQAT(const AccpmGenMatrix &A, const AccpmVector &d, SymmetricMatrix &result) const;
    
    void computeDiagQ(const Manager &manager, const AccpmVector &y, const AccpmGenMatrix *d2f2,
		      double ss, AccpmVector &diagQ) const;
    
    void computeSCBarrier(const Manager &manager, const AccpmVector &y, AccpmVector &result) const;
   
    void computeH(const SymmetricMatrix &H11, const AccpmGenMatrix &H12, const AccpmVector &H13,
		  const SymmetricMatrix &H22, const AccpmVector &H23, double H33,
		  AccpmGenMatrix &H) const;

    void computeFullMatrix(const SymmetricMatrix &H11, AccpmGenMatrix &H) const;

    void updateATQA(Manager &manager, const LocSet &locSet,
		    const AccpmVector &y, const AccpmVector &diagQ) const;
   
    double computeStepSize(const Manager &manager, const LocSet &locSet, 
			   const NewtonSolution &sol, const NewtonSolution &step, 
			   double maxRs, double uNorm) const;
    void computeStepSize(const Manager &manager, const NewtonSolution &sol, const NewtonSolution &step, 
			AccpmVector &S, AccpmVector &DS, AccpmVector &W) const;
    
    double lineSearchDualPot(const Manager &manager, const LocSet &locSet, 
			     const NewtonSolution &sol, const NewtonSolution &step, 
			     const AccpmVector &S, const AccpmVector &DS, const AccpmVector &W) const;

    double maxStep(const AccpmVector &z, const AccpmVector &dz) const;
  
    double maxBall(const Manager &manager, const AccpmVector &y, 
		   const AccpmVector &dy, double alpha) const; //REf: MATLAB: AMaxBall
   
    double computeTerminationCriterion(const Manager &manager, const LocSet &locSet, 
				       const AccpmVector &ry, const AccpmVector &rz,
				       const AccpmVector &wsV, 
				       double ws0, double maxRs, double uNorm, double nnorm, 
				       const AccpmVector *df2, const NewtonSolution &sol, 
				       const NewtonSolution &step) const;
  
    void solveSCSystem(const AccpmVector &a, const AccpmVector &b, const AccpmVector &rhs,
		    AccpmVector &x) const;
    /**
     * Compute Inverse using the Shearman-Morrison Formula
     */
    int computeSCInverse(const AccpmVector &A, const AccpmVector &b, AccpmGenMatrix &Inv) const;
    
    void computeELET(const Manager &manager, const AccpmGenMatrix &E, const AccpmVector &sigma, 
		     AccpmVector &ELET) const;
   
    void handleEqualityConstraints(Manager &manager, LocSet &locSet, 
				   const AccpmVector &wsV, double ws0, 
				   double f2, const AccpmVector *df2,
				   const AccpmGenMatrix *d2f2, 
				   const AccpmVector &diagQ,
				   const NewtonSolution &sol,
				   NewtonSolution &step) const;
    
    int linSolve(Manager &manager, const LocSet &locSet,
			     const AccpmGenMatrix &A, RealMatrix &X,
			     const RealMatrix &B) const;
  public:
    DualMethod(const Parameters *param);
    /**
     * The main loop for the DualMethod.
     * 
     * Each call to this function does the following:
     *
     * 1. Creates a Localization Set, LocSet, based on the 
     * cuurent set of active cuts (and in future variables).
     *
     * 2. Solves the Newton System of equations. It performs Newton Steps,
     * not exceeding the parameter MaxInnerIterations, performNewtonSteps().
     *
     * 3. Computes and update the lower bound, computeLowerBound().
     *
     */
    int run(Manager &manager);

    virtual ~DualMethod();
  };
}
    
#endif
