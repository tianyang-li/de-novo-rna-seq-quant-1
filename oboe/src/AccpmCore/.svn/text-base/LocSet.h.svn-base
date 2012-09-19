// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef LOCSET_H
#define LOCSET_H

/**
 * @file LocSet.h
 * @ingroup AccpmCore
 */
#include "AccpmGenMatrix.h"
#include "Manager.h"

#include "config.h"

#ifdef OBOE_HAS_GLPK
extern "C" {
#include "glpk.h"
}
#endif

#ifdef COIN_HAS_GLPK
extern "C" {
#include "glpk.h"
}
#endif

#ifdef COIN_HAS_OSI
#include "OsiClpSolverInterface.hpp"
#include "OsiGlpkSolverInterface.hpp"
#endif

namespace Accpm 
{

/**
 * Class for the LocalizationSet.
 * 
 * @ingroup AccpmCore */  
  
  class LocSet {

  public:
    LocSet(const Manager &manager, const Parameters &param);
    virtual ~LocSet();
    const AccpmGenMatrix &getA() const { return *_A; }
    const AccpmGenMatrix &getE() const { return _E; }
    const AccpmVector &getC() const { return _c; }
    const AccpmVector &getProximalCenter() const { return _proximalCenter; }
    const double &getRhs() const { return _rhs; }
    const AccpmGenMatrix *getAFull() const { return _AFull; }
    const AccpmGenMatrix *getEFull() const { return _EFull; }
    const AccpmGenMatrix *getATQAFull() const { return _ATQAFull; }
    const AccpmGenMatrix *getAT() const { return _AT; }
    const AccpmGenMatrix *getEFullT() const { return _EFullT; }
    int computeFullAE(const Parameters &param, const AccpmVector *df2);
    int computeFullATQA(const AccpmGenMatrix &ATQA, const AccpmVector &diagQ, const AccpmVector *df2);
    /**
     * Check if the Localization set has a non-empty interior.
     * This solves an LP feasibilty problem using a solver.
     * Currently we use glpk package to solve the LP.
     *
     * @return Returns whether or not the localization set is feasible.
     * If feasible it also returns the simplex objective value in objVal.
     * @returns The primal variables(query point) in y
     * @returns The dual variables in x
     */
    int checkFeasibility(const Parameters &param, double &objVal, AccpmVector &y,
			 AccpmVector &x) const;

  private:
    const AccpmGenMatrix *_A;
    AccpmGenMatrix _E;
    AccpmGenMatrix *_AFull;
    AccpmGenMatrix *_EFull;
    AccpmGenMatrix *_ATQAFull;
    AccpmVector _proximalCenter;
    AccpmVector _c;
    AccpmGenMatrix *_AT;
    AccpmGenMatrix *_EFullT;

    double _rhs;

    void clear();
    /**
     * Helper function for checkFeasibility().
     * Uses Osi interface to solve the LP.
     */
    int solveLP(const Parameters &param, double &objVal, 
		AccpmVector &y, AccpmVector &x) const;
#ifdef COIN_HAS_OSI
    /**
     * Helper function for checkFeasibility().
     * Uses Osi interface with the chosen solver to solve the LP.
     */
    int solveWithOsi(OsiSolverInterface *si, const Parameters &param, 
		      double &objVal,
		      AccpmVector &y, AccpmVector &x) const;
#endif
    /**
     * Helper function for checkFeasibility().
     * Uses glpk package to solve the LP.
     */
    int solveWithGLPK(const Parameters &param, double &objVal, 
		       AccpmVector &y, AccpmVector &x) const;
  };
}

#endif
