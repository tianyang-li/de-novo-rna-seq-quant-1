// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// The OBOE team
//
//

#include "Method.h"
#include "LocSet.h"
#include "Manager.h"
#include "AccpmDefs.h"
#include "Oracle.h"
#include "Solution.h"
#include "AccpmLASolve.h"

#include <cmath>

#include "AccpmBlasInterface.h"


namespace Accpm {

DualMethod::DualMethod(const Parameters *param) : Method(), _param(param)
{
  _computeLB = _param->getIntParameter("ComputeLowerBound");
}
  
int 
DualMethod::run(Manager &manager)
{
  LocSet locSet(manager, *_param);
  manager.getPointGen()._xScaled.resize(1, 0);
  if (performNewtonSteps(manager, locSet) == 0) {
    computeLowerBound(manager, locSet);
  }

  return 0;
}
  
DualMethod::~DualMethod()
{
}

int
DualMethod::performNewtonSteps(Manager &manager, LocSet &locSet)
{
  int my = locSet.getA().size(0);
  int n = locSet.getA().size(1);

  int p = _param->getIntParameter("NumSubProblems");
  NewtonSolution sol(my, n, p);
  sol.init(manager);
  NewtonSolution step(my, n, p);
  
  double nnorm = ACCPM_PLUS_INF;
  _numIter = 0;
  
  AccpmVector sOld;
  double sOld0 = 0;
  double eta = _param->getRealParameter("Eta");
  while (_numIter < _param->getIntParameter("MaxInnerIterations") &&
	 !DBL_LT(nnorm, eta)) {
    ++_numIter;
    
    AccpmVector u = sol._y;
    u.append(sol._z);
    double uNorm = AccpmLANorm2(u);
    AccpmVector rs(n);
    double rs0;
    computeInfeasibilities(manager, locSet, sol._y, sol._z, sol._zs, sol._s, sol._s0, rs, rs0);
    //std::cout << "s =\n" << sol._s << "s0 = " << sol._s0 << std::endl;
    //std::cout << "rs =\n" << rs << "rs0 = " << rs0 << std::endl;
    
    manager.computeSmoothComponent(true);
    double f2 = manager.getSmoothObj();
    if (_param->getIntParameter("Verbosity") > 3) {
      std::cout << "Inner Iteration#:" << _numIter 
		<< " f2: " << f2 << std::endl;
    }
    const AccpmVector *df2 = manager.getSmoothGradient(true);
    const AccpmGenMatrix *d2f2 = manager.getSmoothHessian(true);

    // keep previous iterate and defne fixed values
    double maxRs = AccpmLANormInf(rs);
    if (DBL_LT(maxRs, std::fabs(rs0))) {
      maxRs = std::fabs(rs0);
    }
    sOld = sol._s;
    sOld0 = sol._s0;
    // double ssOld = manager.getPointGen()._ss;
    
    double ws0 = manager.inPhase2() ? (manager.getWtEpigraphCut()*1.0)/sol._s0 : 0;
    AccpmVector wsV = manager.getCutOccurrence();
    wsV.rdivide(sol._s);
    AccpmVector ry(my);
    computeRy(manager, locSet, sol._y, wsV, df2, sol._ss, ry);
    AccpmVector rz(sol._z.size());
    rz = 0;
    if (manager.inPhase2()) {
      computeRz(locSet, wsV, ws0, rz);
    }
    AccpmVector diagQ(my);
    computeDiagQ(manager, sol._y, d2f2, sol._ss, diagQ);
  
    if (newtonDualDir(manager, locSet, ry, rz, rs, rs0, wsV, ws0, f2, df2, d2f2, diagQ, sol, step) == 1) {
      return 1;
    }
   
    if (_param->hasEqualityConstraints()) {
      handleEqualityConstraints(manager, locSet, wsV, ws0, f2, df2, d2f2, diagQ, sol, step);
    }
    
    double alpha = computeStepSize(manager, locSet, sol, step, maxRs, uNorm);
    if (_param->getIntParameter("Verbosity") > 4) {
      std::cout << "Newton Step with alpha:" << alpha << std::endl; step.output(std::cout);
    }
    nnorm = computeTerminationCriterion(manager, locSet, ry, rz, wsV, ws0, maxRs, uNorm, nnorm, df2, sol, step);

    sol.addMult(step, alpha);
    manager.updateY(sol._y);
    if (_param->getIntParameter("Verbosity") > 4) {
      std::cout << "New Solution:\n" << std::endl; sol.output(std::cout);
    }
  }

  manager.updateVariables(sol._y, sol._z, sol._s, sol._zs, sol._s0, sol._ss, 
			  sOld, sOld0, step._s, step._s0);
  if (_param->getIntParameter("Verbosity") > 2) {
    std::cout << "Number of Inner Iterations: " << _numIter << std::endl;
  }
  manager.updateInnerIterations(_numIter);
  return 0;

}

double
DualMethod::computeStepSize(const Manager &manager, const LocSet &locSet,
			    const NewtonSolution &sol, const NewtonSolution &step, 
			    double maxRs, double uNorm) const
{
  AccpmVector S, DS, W;
  computeStepSize(manager, sol, step, S, DS, W);
  double alpha = 1;
  // alpha = min(1, MaxStep(S, DS) * ParamS.Gamma);
  double gamma = _param->getRealParameter("Gamma");
  alpha = std::min(1.0, maxStep(S, DS) * gamma);
  if ((_param->getOracle())->hasSmoothOracle()) {
    if (_param->getIntParameter("Ball")) {
      alpha = maxBall(manager, sol._y, step._y, alpha) * gamma;
    }
  } else {
    if (manager.inPhase2() && DBL_LT(maxRs, _param->getRealParameter("EpsilonReal") * uNorm)) {
      //std::cerr << "Using Line Search Dual Pot" << std::endl;
      alpha = lineSearchDualPot(manager, locSet, sol, step, S, DS, W); 
    } else {
      if (_param->getIntParameter("Ball")) {
	alpha = maxBall(manager, sol._y, step._y, alpha);
      }
    }
  }
  return alpha;
} 

void
DualMethod::computeStepSize(const Manager &manager, const NewtonSolution &sol, const NewtonSolution &step, 
			    AccpmVector &S, AccpmVector &DS, AccpmVector &W) const
{
  /* Bug fix for computeStepSize : July 21, 2008 */
  S.append(sol._s);
  DS.append(step._s);
  W.append(manager.getCutOccurrence());

  if (manager.inPhase2()) {
    S.append(sol._s0);
    S.append(sol._ss);
    DS.append(step._s0);
    DS.append(step._ss);
    W.append(manager.getWtEpigraphCut());
    W.append(manager.getPointGen()._ws);
  }

  if (_param->getIntParameter("Box")) {
    AccpmVector ones(sol._y.size());
    ones = 1;
    const AccpmVector *varB = manager.getVariableLB(true);
    if (varB) {
      AccpmVector tmp = sol._y;
      AccpmLAAddMult(tmp, -1, *varB); 
      S.append(tmp);
      DS.append(step._y);
      W.append(ones);
    }
    varB = manager.getVariableUB(true);
    if (varB) {
      AccpmVector tmp = *varB;
      AccpmLAAddMult(tmp, -1, sol._y); 
      S.append(tmp);
      tmp = step._y;
      tmp.negate();
      DS.append(tmp);
      W.append(ones);
    }
  }
}

/**
 * Compute the residual values for rs and rs0.
 * Update slacks [s,s0] if required.
 */
int 
DualMethod::computeInfeasibilities(const Manager &manager, const LocSet &locSet, 
				   const AccpmVector &y, const AccpmVector &z,
				   double zs,
				   AccpmVector &s, double &s0,
				   AccpmVector &rs, double &rs0)
{
  double tolerance = _param->getRealParameter("EpsilonTol");
  int numCuts = s.size();
  // rs = -(A'y - c + s)
  // rs = -A'*y
  //AccpmLAMatTransVecMult(locSet.getA(), y, rs, -1, 0);
  const AccpmGenMatrix *AT = locSet.getAT();
  AccpmVector ATy(AT->size(0));
  AccpmLAMatVecMult(*AT, y, ATy, 1, 0);
  AccpmLAMult(rs, -1, ATy);
  AccpmLAAddMult(rs, 1, locSet.getC());
  AccpmLAAddMult(rs, -1, s);
  
  if (manager.inPhase2()) {
    AccpmVector Ez(numCuts);
    AccpmLAMatTransVecMult(locSet.getE(), z, Ez);
    // rs = -(A'y - E'z - c + s)
    AccpmLAAddMult(rs, 1, Ez);
    
    // rs0 = - (s0 - LocSetS.rhs + LocSetS.pi' * z + zs);
    double piZ = AccpmLADotProd(*(_param->getPi()), z);
    rs0 = -s0 + locSet.getRhs() - piZ - zs;
    if (DBL_GT(rs0/s0, tolerance)) {
      s0 += rs0;
      rs0 = -s0 + locSet.getRhs() - piZ - zs;
    }
  } else {
    rs0 = 0;
    s0 = 0;
  }
  updateRsS(manager, locSet, ATy, z, rs, s);
  return 0;
}

/**
 * This function updates the values of rs and s depending on the ratio os rs and s
 *
 **/ 
void
DualMethod::updateRsS(const Manager &manager, const LocSet &locSet, 
		      const AccpmVector &ATy, const AccpmVector &z,
		      AccpmVector &rs, AccpmVector &s)
{
  AccpmVector rsbyS = rs;
  rsbyS.rdivide(s);
  bool eEmpty = locSet.getE().size(0) == 0;
 
  const AccpmVector *c = &locSet.getC();
  const AccpmGenMatrix *E = &locSet.getE();
  double tolerance = _param->getRealParameter("EpsilonTol");
  for (int i = 0; i < s.size(); ++i) {
    if (DBL_GT(rsbyS(i), tolerance)) {
      s(i) += rs(i);
      rs(i) = (*c)(i) - ATy(i) - s(i);
      if (!eEmpty) {
	rs(i) +=  AccpmLADotProd(E->getColumn(i), z);
      }
    }
  }
}

int 
DualMethod::newtonDualDir(Manager &manager, LocSet &locSet, 
			  const AccpmVector &ry, const AccpmVector &rz,
			  const AccpmVector &rs, double rs0, 
			  const AccpmVector &wsV, double ws0, 
			  double f2, const AccpmVector *df2,
			  const AccpmGenMatrix *d2f2,
			  const AccpmVector &diagQ,
			  const NewtonSolution &sol,
			  NewtonSolution &step, bool eqConstraint) const
{
  if (!manager.inPhase2()) {
    step._z = 0;
    step._s0 = 0;
  } 
  int my = locSet.getA().size(0);
  int n = locSet.getA().size(1);
  
  double rss = 0;
  double rzs = 0;
  
  if (!eqConstraint) {
    rss = -(sol._ss + f2 - sol._zs);
    rzs = -(ws0 - manager.getPointGen()._ws/sol._ss);
  }
 
  if (n <= _ratioCutVsDim * my) {
    //std::cerr << "Newton Dual Direction for case: NumCuts <= NumVariables" << std::endl;
    return newtonDualDir1(manager, locSet, ry, rs, rs0, rz, rzs, rss, diagQ, df2, sol, step, eqConstraint);
    
  } else { /* Number of cuts > ratio * problem dimension */
    //std::cerr << "Newton Dual Direction for case: NumCuts > NumVariables" << std::endl;
    return newtonDualDir2(manager, locSet, ry, rs, rs0, wsV, rz, rzs, rss, diagQ, df2, sol, step);
  }

  return 0;
}

/**
 * Newton Step for the case Num of Cuts > Number of Variables
 *
 */
int 
DualMethod::newtonDualDir2(Manager &manager, LocSet &locSet, 
			   const AccpmVector &ry, 
			   const AccpmVector &rs, double rs0,
			   const AccpmVector &wsV,
			   const AccpmVector &rz, double rzs, double rss,
			   const AccpmVector &diagQ,
			   const AccpmVector *df2,
			   const NewtonSolution &sol,
			   NewtonSolution &step) const
{
  int status = 0;
  int my = locSet.getA().size(0);
  
  AccpmVector sigma = wsV;
  sigma.rdivide(sol._s); 
  SymmetricMatrix H(my, my);
  prodAQAT(locSet.getA(), sigma, H);
  for (int i = 0; i < my; H(i,i) += diagQ(i), ++i);
  if (_param->getIntParameter("Ball")) {
    AccpmVector rang(my);
    manager.computeBallConstraint(sol._y, rang);
    AccpmLAR1Update(H, rang, 1);
  }
  //rhs = ry + LocSetS.A * Sigma * rs;
  AccpmGenMatrix Ad(locSet.getA());
  for (int i = 0; i < Ad.size(1); ++i) {
    Ad.scaleColumn(i, sigma(i));
  }
  AccpmVector rhs = ry;
  AccpmLAMatVecMult(Ad, rs, rhs, 1, 1);
  
  if (manager.inPhase2()) {
    double w0 = manager.getWtEpigraphCut();
    double w0Bys0Sq = w0 / pow(sol._s0, 2);
    double ws = manager.getPointGen()._ws;
    double wsByssSq = ws / pow(sol._ss, 2);
    AccpmVector rhs1 = rhs;
    if (df2) {
      AccpmLAR1Update(H, *df2, wsByssSq);
      // rhs1 = rhs + ws * rss * df2 / ss ^ 2;
      AccpmLAAddMult(rhs1, wsByssSq*rss, *df2);
    }
    // rhs2 = rz - LocSetS.E * Sigma * rs + LocSetS.pi * rs0 * w0 / s0.^2;
    AccpmVector rhs2 = rz;
    AccpmGenMatrix Ed(locSet.getE());
    for (int i = 0; i < Ed.size(1); ++i) {
      Ed.scaleColumn(i, sigma(i));
    }
    AccpmLAMatVecMult(Ed, rs, rhs2, -1, 1);
    double scale = rs0 * w0Bys0Sq;
    AccpmLAAddMult(rhs2, scale, *_param->getPi());
    
    // rhs3 = rzs + rs0 * w0 / s0.^2 - ws * rss / ss ^ 2;
    double rhs3 = rzs + rs0 * w0Bys0Sq - rss * wsByssSq;
   
    int p = locSet.getE().size(0);
    // H12 = - LocSetS.A * Sigma * LocSetS.ET;
    /* Compute the transpose as we like to keep the lower triangular matrix */
    AccpmGenMatrix H21(p, my);
    AccpmLAMatMatTransMult(locSet.getE(), Ad, H21, -1, 0);
   
    // H13 = - ws * df2 / ss ^ 2;
    AccpmVector H13(my);
    if (df2) {
      AccpmLAMult(H13, -wsByssSq, *df2); 
    } else {
      H13 = 0;
    }
    // H22 = LocSetS.E * Sigma * LocSetS.ET + LocSetS.pi * LocSetS.pi' * w0 / s0.^2;
    SymmetricMatrix H22(p, p);
    H22 = 0;
    AccpmVector ELET(p);
    computeELET(manager, locSet.getE(), sigma, ELET);
    
    if (p > 1) {
      status = newtonDualDirForLargeP(manager, locSet,
				      H, H21, H13, ELET, wsByssSq, w0Bys0Sq, rhs1, rhs2, rhs3, step);
    } else {
      for (int i = 0; i < p; ++i) {
	H22(i,i) = ELET(i);
      }
      AccpmLAR1Update(H22,  *_param->getPi(), w0Bys0Sq);
      
      // H23 = LocSetS.pi * w0 / s0.^2;
      AccpmVector H23(p);
      AccpmLAMult(H23, w0Bys0Sq, *_param->getPi());
      
      // H33 = ws / ss^2 + w0 / s0.^2;
      double H33 = wsByssSq + w0Bys0Sq;
      
      // H = [ H H12 H13; H12' H22 H23; H13' H23' H33];
      AccpmGenMatrix CH(my + p + 1, my + p + 1);
      computeH(H, H21, H13, H22, H23, H33, CH);
      AccpmVector v(my+p+1);
   
      AccpmVector rhs4;
      rhs4.append(rhs1);
      rhs4.append(rhs2);
      rhs4.append(rhs3);
      //std::cout << "rhs4: " << rhs4 << std::endl;
      
      status = linSolve(manager, locSet, CH, v, rhs4);
      if (status == 0) {
	step._y = v(LaIndex(0,my-1));
	
	//  step._z = v(LaIndex(my, my+p-1)); // does not work because of new copy
	for (int i = 0; i <= p-1; ++i) {
	  step._z(i) = v(i+my);
	}
	step._zs = v(my+p);
	
      } else {
	return status;
      }
    }
    if (status == 0) {
      // ds = rs - LocSetS.AT * dy;
      step._s = rs;
      //AccpmLAMatTransVecMult(locSet.getA(), step._y, step._s, -1, 1);
      AccpmLAMatVecMult(*locSet.getAT(), step._y, step._s, -1, 1);
      step._s0 = rs0 - step._zs - AccpmLADotProd(*_param->getPi(), step._z);
      //ds = ds + FastProd(LocSetS.ET, dz, LocSetS.E);
      AccpmLAMatTransVecMult(locSet.getE(), step._z, step._s, 1, 1);
      //dss = rss - df2' * dy + dzs;
      step._ss = rss + step._zs;
      if (df2) {
	step._ss -= AccpmLADotProd(*df2, step._y);
      }
    }
  } else {
    AccpmVector v(my);
    if (AccpmLASymmLinSolve(H, v, rhs) != 0) {
      manager.setExitCode(LA_ERROR);
    }
    
    step._y = v(LaIndex(0,my-1));
    // ds = rs - LocSetS.AT * dy;
    step._s = rs;
    //AccpmLAMatTransVecMult(locSet.getA(), step._y, step._s, -1, 1);
    AccpmLAMatVecMult(*locSet.getAT(), step._y, step._s, -1, 1);
  }
  
  return status;
}

/**
 * Computes EWS-2ET matrix which is a diagonal matrix since E is orthogonal.
 **/
void
DualMethod::computeELET(const Manager &manager, const AccpmGenMatrix &E, const AccpmVector &sigma, 
			AccpmVector &ELET) const
{
  ELET = 0;
  const AccpmVector &subProblemIndex = manager.getSubProblemIndex();
  assert(subProblemIndex.size() == sigma.size());
  int numCuts = subProblemIndex.size();
  for (int i = 0; i < numCuts; ++i) {
    int problemIndex = (int)subProblemIndex(i);
    if (problemIndex > 0) {
      ELET(problemIndex - 1) += sigma(i);
    }
  }
}
/**
 * Newton Step for the case Num of Cuts <= Number of Variables
 *
 */
int 
DualMethod::newtonDualDir1(Manager &manager, LocSet &locSet, 
			   const AccpmVector &ry, 
			   const AccpmVector &rs, double rs0,
			   const AccpmVector &rz, double rzs, double rss,
			   const AccpmVector &diagQ,
			   const AccpmVector *df2,
			   const NewtonSolution &sol,
			   NewtonSolution &step,
			   bool eqConstraint) const
{
  int my = locSet.getA().size(0);
  int n = locSet.getA().size(1);
  AccpmVector sigma = manager.getCutOccurrence();
  AccpmVector sigmaMinus;
  AccpmGenMatrix H;
  AccpmVector Qry(ry);
  Qry.rdivide(diagQ);

  AccpmVector ATQry(n);

  AccpmVector rs1(rs);
     
  AccpmVector rz1(rz);
  if (manager.inPhase2()) {
    rs1.append(rs0);
    rs1.append(rss);
    rz1.append(rzs);
  }

  if (!eqConstraint) {
    if (_numIter == 1 || _param->getIntParameter("Ball") ||
	_param->getIntParameter("Box") || 
	(_param->getOracle())->hasSmoothOracle()) {
      updateATQA(manager, locSet, sol._y, diagQ);
    }
  }
  if (manager.inPhase2()) {
    if (!eqConstraint) {
      if (_numIter == 1 ||  (_param->getOracle())->hasSmoothOracle()) {
	// LocSetS.AFull = [LocSetS.A zeros(my,1) df2];
	locSet.computeFullAE(*_param, df2);
	// LocSetS.EFull = [LocSetS.E  -LocSetS.pi zeros(mz,1); zeros(1,n) -1 1];
      }
      if (_numIter == 1 || _param->getIntParameter("Ball") ||
	  _param->getIntParameter("Box") || 
	  (_param->getOracle())->hasSmoothOracle()) {
	locSet.computeFullATQA(manager.getATQA(), diagQ, df2);
      }
    }
   
    double w0 = manager.getWtEpigraphCut();
    double ws = manager.getPointGen()._ws;
    sigma.append(w0);
    sigma.append(ws);
    AccpmVector s; 
    s.append(sol._s);
    s.append(sol._s0);
    s.append(sol._ss);
    s.times(s);
    sigma.rdivide(s);
    sigmaMinus = sigma;
    sigmaMinus.invert();
    H = *locSet.getATQAFull();
    ATQry.resize(n+2, 1);
    AccpmLAMatTransVecMult(*locSet.getAFull(), Qry, ATQry, 1, 0);
  } else {
    sigma.rdivide(sol._s);
    sigma.rdivide(sol._s);
    sigmaMinus = sigma;
    sigmaMinus.invert();
    H = manager.getATQA();
    //AccpmLAMatTransVecMult(locSet.getA(), Qry, ATQry, 1, 0);
    AccpmLAMatVecMult(*locSet.getAT(), Qry, ATQry, 1, 0);
  }

  for (int i = 0; i < H.size(0); H(i,i) += sigmaMinus(i), ++i);
  AccpmVector rhs = rs1;
  
  AccpmLAAddMult(rhs, -1, ATQry);
  
  AccpmGenMatrix L;
  bool cholesky = true;
  if (AccpmLACholeskyFactor(H, L) > 0) {
    manager.setExitCode(CHOLESKY_FAILURE);
    if (_param->getIntParameter("CheckLocSetInterior")) {
      double kellyBound;
      AccpmVector ky;
      AccpmVector kx;
      int feasible = locSet.checkFeasibility(*_param, kellyBound, ky, kx);
      if (feasible != 1) {
	manager.setExitCode(LOCSET_EMPTY);
	return 1;
      } else if (feasible == 1) {
	if (IS_FINITE(kellyBound)) {
	  std::cout << "Updating lower bound with the obtained Kelly bound: " 
		    << kellyBound << std::endl;
	  manager.updateLB(kellyBound);
	  manager.getPointGen()._xScaled = kx;
	  manager.setExitCode(ITERATING);
	  return 1;
	}
      }
    }
    cholesky = false;
    L.copy(H);
  }
  if (manager.inPhase2()) {
    AccpmVector nrhs(rhs);
    nrhs.negate();
    AccpmVector tmpresult(H.size(1));
    // dz = rz + E* LinSolve(H, R, -rhs);
    //LaLinearSolve(H, tmpresult, nrhs);
    
    AccpmLALinSolve(L, cholesky, tmpresult, nrhs);
    step._z = rz1;
    AccpmLAMatVecMult(*locSet.getEFull(), tmpresult, step._z, 1, 1);
    
    // dz = (E * LinSolve(H, R,ET)) \ dz;
    const AccpmGenMatrix *E = locSet.getEFull();
    const AccpmGenMatrix *ET = locSet.getEFullT();
     
    int mp = ET->size(1);
    assert(mp == 1+_param->getIntParameter("NumSubProblems"));

    AccpmGenMatrix tmp1(H.size(1), mp);
    AccpmLALinSolve(L, cholesky, tmp1, *ET);
    
    AccpmGenMatrix tmp2(mp, mp);
    AccpmLAMatMatMult(*E, tmp1, tmp2, 1, 0);
    tmpresult.resize(mp, 1);
    AccpmLALinearSolve(tmp2, tmpresult, step._z);
    step._z = tmpresult;
     
   // rhs = rhs + ET * dz;
    AccpmLAMatVecMult(*ET, step._z, rhs, 1, 1);

    step._zs = step._z(mp-1);
    //step._z = step._z(LaIndex(0, mp-2));
    AccpmVector ztmp(mp-1);
    for(int i = 0; i < mp-1; ++i) {
      ztmp(i) = step._z(i);
    }
    step._z = ztmp;
  }
  // ds = SigmaMinus * LinSolve(H, R,rhs);
  AccpmVector result(H.size(1));
 
  AccpmLALinSolve(L, cholesky, result, rhs);
  step._s = sigmaMinus;
  step._s.times(result);

   //dy  = diaQ .\ (FastProd(A, (Sigma * ds), AT) + ry);
  sigma.times(step._s);
  AccpmVector result1 = ry;

  if (manager.inPhase2()) {
    AccpmLAMatVecMult(*locSet.getAFull(), sigma, result1, 1, 1);
    
  } else {
    AccpmLAMatVecMult(locSet.getA(), sigma, result1, 1, 1);
  }
  //std::cout << "Result1:\n" << result1 << std::endl;
  if (_param->getIntParameter("Ball")) {
    // SolveSCSystem?
    AccpmVector B(my);
    manager.computeBallConstraint(sol._y, B);
    AccpmVector tmp = diagQ;
    tmp.negate();
    
    solveSCSystem(tmp, B, result1, step._y);
    step._y.negate();
  
  } else {
    step._y = result1;
    step._y.rdivide(diagQ);
  }

  if (manager.inPhase2()) {
    step._ss = step._s(n+1);
    step._s0 = step._s(n);
    step._s = step._s(LaIndex(0, n-1));
  }

  return 0;
}

int
DualMethod::computeRy(const Manager &manager, const LocSet &locSet, const AccpmVector &y, const AccpmVector &wsV, 
		      const AccpmVector *df2, double ss,
		      AccpmVector &ry) const
{
  // ry = -(Q(y-y') + wf2'(y/)sigma +H'(y) +AWs-1)
  AccpmVector tmp = y;
  AccpmLAAddMult(tmp, -1, locSet.getProximalCenter());
  tmp.times(*manager.getQ(true));
  // Q(y-y') + Aws-1
  AccpmLAMatVecMult(locSet.getA(), wsV, tmp, 1, 1);
  // + H'(y)
  
  manager.computeBallConstraint(y, ry);
  AccpmLAAddMult(ry, 1, tmp);
  
  manager.computeBox1Constraint(y, tmp);
  AccpmLAAddMult(ry, 1, tmp);

  // +wf2'(y)/sigma  
  if (manager.inPhase2() && df2) {
    double scale = manager.getPointGen()._ws/ss;
    AccpmLAAddMult(ry, scale, *df2); 
  }

  ry.negate();
  
  return 0;
}

int
DualMethod::newtonDualDirForLargeP(Manager &manager, LocSet &locSet,
				   const SymmetricMatrix &H, const AccpmGenMatrix &H21, 
				   const AccpmVector &H13, const AccpmVector &ELET, double wsByssSq, 
				   double w0Bys0Sq, const AccpmVector &rhs1, const AccpmVector &rhs2, double rhs3, 
				   NewtonSolution &step) const
{
  int status = 0;
  AccpmVector A(ELET);
  A.append(wsByssSq);
  AccpmVector b(*(_param->getPi()));
  b.append(1);
  AccpmLAScale(sqrt(w0Bys0Sq), b);
  int p = b.size(); // Number of subproblems + 1
  AccpmGenMatrix LInv(p, p);
  computeSCInverse(A, b, LInv);
 
  int my = H.size(0);
  AccpmGenMatrix B(my, p); //B = [H12 H13] = [H21^t  H13]
  for (int i = 0; i < p-1; ++i) {
    B.assignColumn(i, H21.getRow(i));
  }
  B.assignColumn(p-1, H13);
  
  AccpmGenMatrix BLInv(my, p);
  AccpmLAMatMatMult(B, LInv, BLInv, 1, 0);
 
  AccpmGenMatrix Y(my, my);
  computeFullMatrix(H, Y);
  AccpmLAMatMatTransMult(BLInv, B, Y, -1, 1);

  AccpmVector rhs;
  rhs.append(rhs2);
  rhs.append(rhs3);
  AccpmVector result(rhs1);
  //rhs1 - B*LInv*rhs
  AccpmLAMatVecMult(BLInv, rhs, result, -1, 1);
  
  status = linSolve(manager, locSet, Y, step._y, result);
  if (status == 0) {
    //dz = LInv * rhs - LInv*B^t*dy
    AccpmVector LInvrhs(p);
    AccpmLAMatVecMult(LInv, rhs, LInvrhs, 1, 0);
    AccpmVector dz(p);
    AccpmLAMatTransVecMult(BLInv, step._y, dz, -1, 0);
    AccpmLAAddMult(dz, 1, LInvrhs);
    step._z = dz(LaIndex(0, p-2));
    step._zs = dz(p-1);
  }

  return status;
}

/**
 * This function computes the inverser for matrix A + bb^t
 * via Sherman-Morrison formula. The matrix A is in the form of non-zero diagonal vector.
 */
int
DualMethod::computeSCInverse(const AccpmVector &A, const AccpmVector &b, AccpmGenMatrix &Inv) const
{
  AccpmVector bbyA(b);
  bbyA.rdivide(A);
  AccpmVector nbbyA(bbyA);
  nbbyA.negate();
  double factor = 1.0/(1 + AccpmLADotProd(bbyA, b));
  AccpmLAScale(factor, bbyA);
  
  Inv = 0;
  AccpmVector AInv(A);
  AInv.invert();
  for (int i = 0; i < A.size();  Inv(i,i) += AInv(i), ++i);

  AccpmLAR1Update(Inv, bbyA, nbbyA, 1);

#if 0
  std::cout << "A:\n" << A << "b:\n" << b << std::endl;
  std::cout << "Inverse for [A bbT]:\n" << Inv << std::endl; 
#endif

  return 0;
}

int
DualMethod::computeRz(const LocSet &locSet, const AccpmVector &wsV, double ws0, AccpmVector &rz) const
{
  // rz = -(-EWs-1 + pi*w0s0-1
  // -pi*w0s0-1
  rz = *(_param->getPi());
  AccpmLAScale(-ws0, rz);
  // + EWs-1
  AccpmLAMatVecMult(locSet.getE(), wsV, rz, 1, 1);
 
  //std::cout << "rz:\n" << rz << std::endl;
  return 0;
}

/**
 * compute AdAt and store the result in Symmetric Matrix result.
 */
int
DualMethod::prodAQAT(const AccpmGenMatrix &A, const AccpmVector &d, SymmetricMatrix &result) const
{
  // Scale the columns of A by sqrt d
  assert(A.size(1) == d.size());
  AccpmGenMatrix Ad = A;
  for (int i = 0; i < A.size(1); ++i) {
    Ad.scaleColumn(i, sqrt(d(i)));
  }
  AccpmLAR1Update(result, Ad, 1, 0);
#if 0
  std::cout << "A:\n" << A << std::endl;
  std::cout << "d:\n" << d << std::endl;
  std::cout << "Ad:\n" << Ad << std::endl;
  std::cout << "ProdAQAT:\n" << result << std::endl;
#endif
 return 0;    
}
/**
 * compute the Proximal matrix: Q += H'(y) + ws*ss-1*d2f2
 **/
void
DualMethod::computeDiagQ(const Manager &manager, const AccpmVector &y, const AccpmGenMatrix *d2f2,
			 double ss, AccpmVector &diagQ) const
{
  diagQ = *manager.getQ(true);
  AccpmVector result(diagQ.size());
  computeSCBarrier(manager, y, result);
  AccpmLAAddMult(diagQ, 1, result);
  if (d2f2) {
    if (_param->getIntParameter("DiagHessian")) {
      double factor = manager.getPointGen()._ws/ss;
      AccpmLAAddMult(diagQ, factor, *d2f2);
    } else {
      std::cout << "General Hessian not yet supported: Only Diagonal Hessian allowed" << std::endl;
      exit(0);
    }
  }
}

/**
 * Compute the contribution of the Self-Concordant Barrier H'(y) function
 **/
void
DualMethod::computeSCBarrier(const Manager &manager, const AccpmVector &y, AccpmVector &result) const
{ 
  result = 0;
  if (_param->getIntParameter("Box")) {
    manager.computeBox2Constraint(y, result);
  }
  if (_param->getIntParameter("Ball")) {
    double bc = manager.computeBallConstraint(y);
    AccpmVector tmp(result.size());
    tmp = bc;
    AccpmLAAddMult(result, 1, tmp);
  }
}

/**
 * Compute the complete H matrix from the components
 * Only the lower triangular part of the symmetric matrix is computed
 * and stored in AccpmGenMatrix H.
 */
void
DualMethod::computeH(const SymmetricMatrix &H11, const AccpmGenMatrix &H21, const AccpmVector &H13,
		     const SymmetricMatrix &H22, const AccpmVector &H23, double H33,
		     AccpmGenMatrix &H) const 
{
  int my = H11.size(0);
  int p = H21.size(0);
  LaIndex index1 = H11.index(0);
  //  H(index1, index1) = H11;
  for (int i = 0; i < my; ++i) {
    for (int j = 0; j <= i; ++j) {
      H(i,j) = H11(i,j);
      //H(j,i) = H11(i,j);
    }
  }
  LaIndex rowindex2(my, my+p-1);
  H(rowindex2, index1).inject(H21);
  

  //H(rowindex2, rowindex2) = H22; 
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j <= i; ++j) {
      H(my+i,my+j) = H22(i,j);
      //H(my+j,my+i) = H22(i,j);
    }
  }
  AccpmVector lastRow;
  lastRow.append(H13);
  lastRow.append(H23);
  lastRow.append(H33);
  H.assignRow(my+p, lastRow);
  //H.assignColumn(my+p, lastRow);
#if 0
  std::cout << "H11:  " << H11 << std::endl;
  std::cout << "H21:  " << H21 << std::endl;
  std::cout << "H13:\n" << H13 << std::endl;
  std::cout << "     H22:\n" << H22 << std::endl;
  std::cout << "            H23:\n" << H23 << std::endl;
  std::cout << "                H33:\n\n" << H33 << std::endl;
  std::cout << "H:\n" << H << std::endl;
#endif
}
/**
 * Helper function to transform a SymmetricMatrix into AccpmGenMatrix.
 * This is required since LAPACK++ currently does not support
 * solving AX = B for Symmetric Matrices.
 *
 */
void
DualMethod::computeFullMatrix(const SymmetricMatrix &H11, AccpmGenMatrix &H) const
{
  int my = H11.size(0);
  for (int i = 0; i < my; ++i) {
    for (int j = 0; j <= i; ++j) {
      H(i,j) = H11(i,j);
      //H(j,i) = H11(i,j);// BLAS routines only reference the lower triangular matrix
    }
  }
}

void 
DualMethod::updateATQA(Manager &manager, const LocSet &locSet,
		       const AccpmVector &y, const AccpmVector &diagQ) const
{
  int n = locSet.getA().size(1);
  int sizeOfs = manager.getSizeOfs();
  if (sizeOfs == 0 || manager.getATQA().size(0) == 0 
      || _param->getIntParameter("Ball") ||
      _param->getIntParameter("Box") || 
      (_param->getOracle())->hasSmoothOracle()) {
    manager.getPointGen().setRho(manager.getRho());
    AccpmVector invQ = diagQ;
    invQ.invert();
    int n = locSet.getA().size(1);
    SymmetricMatrix H(n, n);
    const AccpmGenMatrix *AT = locSet.getAT();
    prodAQAT(*AT, invQ, H);
    if (manager.getATQA().size(0) !=  n) {
      manager.getATQA().resize(n, n);
    }
    computeFullMatrix(H, manager.getATQA());
  } else {
    
    AccpmGenMatrix *ATQA = &(manager.getATQA());
    //%If rho is dynamic we rescale ATQA.
    //%ATQA = (RhoOld / Rho)  * ATQAold.
    //if ParamS.DynamicRho == 1
    //aux = PointGenS.RhoOld / ManagerS.D.Rho;
    //PointGenS.ATQA = aux * PointGenS.ATQA;
    //PointGenS.RhoOld = ManagerS.D.Rho;
    //end
    if (_param->getIntParameter("DynamicRho")) {
      double aux = manager.getPointGen()._rho / manager.getRho();
      manager.getATQA().scale(aux);
      manager.getPointGen()._rho =  manager.getRho();
    }
    int numNewCuts = n - sizeOfs;
    assert(numNewCuts >= 0);
    if (numNewCuts) {
      //A1  = LocSetS.A(: , 1:ManagerS.D.SizeOfs);
      LaIndex rowIndex = locSet.getA().index(0);
      LaIndex colIndex(0, sizeOfs-1);
      AccpmGenMatrix A1(locSet.getA()(rowIndex, colIndex));
      //B   = LocSetS.A(: , ManagerS.D.SizeOfs + 1 : n);
      colIndex = LaIndex(sizeOfs, n-1);
      AccpmGenMatrix B(locSet.getA()(rowIndex, colIndex));
      //BTQ = B' * spdiags(ones(my,1) ./ PointGenS.q, 0, my, my);
      AccpmGenMatrix *BTQ = B.transpose();
      for (int i = 0; i < BTQ->size(1); ++i) {
	BTQ->scaleColumn(i, 1/diagQ(i));
      }
      
      //% Compute ATQA = [A1TQA1 A1TQB ; BTQA1 BTQB]
      AccpmGenMatrix BTQB(numNewCuts, numNewCuts);
      AccpmLAMatMatMult(*BTQ, B, BTQB, 1, 0);
      //BTQA = BTQ * A1 ;
      AccpmGenMatrix BTQA(numNewCuts, sizeOfs);
      AccpmLAMatMatMult(*BTQ, A1, BTQA, 1, 0);
      
      //PointGenS.ATQA = Concat( ...
      //PointGenS.ATQA(1:ManagerS.D.SizeOfs, ...
      //1 : ManagerS.D.SizeOfs), BTQA', BTQA, BTQB );
      AccpmGenMatrix newATQA(n, n);
      newATQA(ATQA->index(0), ATQA->index(1)).inject(*ATQA);
      
      for(int i = 0; i < numNewCuts; ++i) {
	AccpmVector v;
	v.append(BTQA.getRow(i));
	v.append(BTQB.getRow(i));
	newATQA.assignRow(sizeOfs+i, v);
      }
      //ATQA->ref(newATQA);
      manager.getATQA() = newATQA;
      delete BTQ;
    }
  }
  manager.setSizeOfs(locSet.getA().size(1));
}

double
DualMethod::maxStep(const AccpmVector &z, const AccpmVector &dz) const
{
  // MaxStep (z,dz) computes the largest amax such that z+amax*dz >= 0.
  assert(z.size() == dz.size());
  double amax = ACCPM_PLUS_INF;
  for (int i = 0; i < dz.size(); ++i) {
    if (DBL_LT(dz(i), 0)) {
      double val = -z(i)/dz(i);
      if (DBL_LT(val, amax)) {
	amax = val;
      }
    } 
  }
  return amax;
}

double
DualMethod::maxBall(const Manager &manager,const AccpmVector &y, 
		    const AccpmVector &dy, double alpha) const
{
  double returnValue = alpha;
  // Computes the largest amax for the ball
  double radius2 = pow(manager.getRadiusBall(), 2);
  AccpmVector centeredY = y;
  AccpmLAAddMult(centeredY, -1, *manager.getCenterBall(true));
  double value1 = radius2 - AccpmLADotProd(centeredY, centeredY);
  double dyTcy = AccpmLADotProd(dy, centeredY);
  double dyTdy = AccpmLADotProd(dy, dy);
  double value2 = pow(alpha, 2) * dyTdy + 2 * alpha * dyTcy;
  if (DBL_LT(value1, value2)) {
    double b = 2 * dyTcy;
    double det = pow(b, 2) + 4 * dyTdy * value1;
    returnValue = (- b + sqrt(det)) / (2 * dyTdy );
  }
  return returnValue;
}
 
double
DualMethod::computeTerminationCriterion(const Manager &manager, const LocSet &locSet,
					const AccpmVector &ry, const AccpmVector &rz,
					const AccpmVector &wsV, 
					double ws0, double maxRs, double uNorm, double nnorm, 
					const AccpmVector *df2, const NewtonSolution &sol, 
					const NewtonSolution &step) const
{
  double returnValue = nnorm;
  if (DBL_LT(maxRs, _param->getRealParameter("EpsilonReal") * uNorm)) {
    if (manager.inPhase2()) {
      if ((_param->getOracle())->hasSmoothOracle()) {
	double rzs = -(ws0 - manager.getPointGen()._ws/sol._ss);
	returnValue = fabs(AccpmLADotProd(ry, step._y) + AccpmLADotProd(rz, step._z) + rzs*step._zs);
      } else {
	returnValue = fabs(AccpmLADotProd(ry, step._y) + AccpmLADotProd(rz, step._z));
      }
    } else {
      returnValue = fabs(AccpmLADotProd(ry, step._y));
    }
  }
  return returnValue;
}

/* Performs a line search on the dual proximal barrier
** DualProxBar = 0.5 a \alpha^2  + b \alpha
* - sum (w.*log(1 + \alpha c)) - log (R2 - norm(y - yc)^2).
*/

double
DualMethod::lineSearchDualPot(const Manager &manager, const LocSet &locSet, 
			      const NewtonSolution &sol, const NewtonSolution &step, 
			      const AccpmVector &sFull, const AccpmVector &dsFull, const AccpmVector &w) const
{
  // a = dy' * (PointGenS.q .* dy);
  AccpmVector qdy = step._y;
  qdy.times(*manager.getQ(true));
  double a = AccpmLADotProd(step._y, qdy);
  // b = dy' * (PointGenS.q .* (y - LocSetS.ProxCenter(1:my)));
  AccpmVector centeredY = sol._y;
  AccpmLAAddMult(centeredY, -1, locSet.getProximalCenter());
  centeredY.times(*manager.getQ(true));
  double b = AccpmLADotProd(step._y, centeredY);
  
  // c = DSfull ./ Sfull;
  AccpmVector c = dsFull;
  c.rdivide(sFull);
  double amax = maxStep(sFull, dsFull);
   
  if (_param->getIntParameter("Ball")) {
    amax = maxBall(manager, sol._y, step._y, amax);
  } 
  int maxIter = 100;
  int nbIter = 0;
  double lNorm = ACCPM_PLUS_INF;
  double precision = 1e-8;
  double alpha = std::min(1.0, 0.9 * amax);
 
  if (_param->getIntParameter("Ball")) {
    centeredY = sol._y;
    AccpmLAAddMult(centeredY, -1, *manager.getCenterBall(true));
  }
  double dySq = pow(AccpmLANorm2(step._y), 2);
  double alphaOld = alpha;
  
  while (DBL_GT(lNorm, precision) && (nbIter < maxIter)) {
    alphaOld = alpha;
    //d = c ./ (1 + alpha * c);
    AccpmVector d = c;
    AccpmVector tmp(c.size());
    tmp = 1;
    AccpmLAAddMult(tmp, alpha, c);
    d.rdivide(tmp);
    //dDualProxBar = a * alpha + b - sum(w .* d);
    double dDualProxBar = a *alpha + b - AccpmLADotProd(w, d); 
    //d2DualProxBar = a + sum(w .* d .* d);
    d.times(d); 
    double d2DualProxBar = a + AccpmLADotProd(w, d);
    if (_param->getIntParameter("Ball")) {
      // g = y + alpha * dy - ParamS.CenterBall;
      AccpmVector g = centeredY;
      AccpmLAAddMult(g, alpha, step._y);
      // e = 2 * dy' * g;
      double e = 2 * AccpmLADotProd(step._y, g);
      // f = ParamS.RadiusBall ^ 2 - norm(g) ^ 2;
      double f = manager.getRadiusBall();
      f = pow(f,2) - pow(AccpmLANorm2(g), 2);
      dDualProxBar +=  e / f;
      //d2DualProxBar = d2DualProxBar + 2 * dy' * dy / f + e * e' / f ^ 2;
      d2DualProxBar += 2*dySq/f + pow(e, 2)/pow(f,2);
    }
    //LNorm = sqrt(d2DualProxBar \ dDualProxBar ^ 2);
    lNorm = std::fabs(dDualProxBar) / sqrt(d2DualProxBar);

    // Damped or Full Newton step (minimizing a self-concordant function)
    if (DBL_GT(3*lNorm, 1.0)) {
      alpha = alpha - (dDualProxBar / d2DualProxBar) / (1 + lNorm);
    } else {
      alpha = alpha - (dDualProxBar / d2DualProxBar);
    }
    ++nbIter;
  }
  // get alpha for which termination criteria was evaluated
  alpha = alphaOld;
  if (nbIter == maxIter) {
    AccpmWarning("LineSearchDualPot : maximum number of iterations reached");
  }
  return alpha;
}

double
DualMethod::computeLowerBound(Manager &manager, const LocSet &locSet)
{
  if (!manager.inPhase2()) {
    manager.getPointGen()._xScaled = 0;
    return manager.getObjLB();
  }
  int my = locSet.getA().size(0);
  int n = locSet.getA().size(1);
  int p = locSet.getE().size(0);
  double lBound = ACCPM_MINUS_INF;

  if (IS_FINITE(manager.getObjUB())) {
    AccpmVector xOld = manager.getCutOccurrence();
    xOld.append(manager.getWtEpigraphCut());
    xOld.rdivide(manager.getPointGen()._sOld);
    //x = xOld .* (1 - PointGenS.DecVar.ds ./ PointGenS.DecVar.sOld);
    AccpmVector tmp = manager.getPointGen()._ds;
    tmp.rdivide(manager.getPointGen()._sOld);
    AccpmVector x(tmp.size());
    x = 1;
    AccpmLAAddMult(x, -1, tmp);
    x.times(xOld);

    if (DBL_LT(x.min(), 0)) {
      x = manager.getPointGen()._x;
    }
    // scaling factor to build a (Wolfe) primal feasible x for original problem.
    double scaleFactor = 1.0/x(n);
    AccpmVector xScaled;
    xScaled = x(LaIndex(0, n-1));
    AccpmLAScale(scaleFactor, xScaled);
    
    if (_computeLB) {
      const AccpmVector &y = manager.getActiveY(); 
      AccpmVector yDelta = y;
      AccpmLAAddMult(yDelta, -1, locSet.getProximalCenter());
      AccpmVector scaledQyDelta = *manager.getQ(true);
      scaledQyDelta.times(yDelta);
      AccpmLAScale(scaleFactor, scaledQyDelta);
      double startDualLB = AccpmLANorm2(scaledQyDelta) / std::max(fabs(manager.getObjUB()), 1.0);
      if (DBL_LT(startDualLB, _startLB)) {
	AccpmVector ExScaled(p);
	AccpmLAMatVecMult(locSet.getE(), xScaled, ExScaled, 1, 0);
        AccpmVector scale = *_param->getPi();
	scale.rdivide(ExScaled);

        // Scale the convexity rows corresponding to the subproblems
        // xScaledTemp = sum(sparse(1:mz,1:mz,scale) * (LocSetS.E .* repmat(xScaled',mz,1)),1);
	AccpmVector xScaledTemp(n);
	AccpmVector v;
	AccpmGenMatrix *ET = locSet.getE().transpose(); // We cant scale rows so take the transpose
	AccpmVector ones(p);
	ones = 1;
	for (int i = 0; i < p; ++i) {
	  ET->scaleColumn(i, xScaled);
	  ET->scaleColumn(i, scale(i));
	}

	for (int i = 0; i < n; ++i) {
	  AccpmVector tmp = ET->getRow(i);
	  xScaledTemp(i) = AccpmLADotProd(tmp, ones);
	}

	delete ET;

	// posindex = find(xScaledTemp > 0);
	// xScaled(posindex) = xScaledTemp(posindex);
	for (int i = 0; i < n; ++i) {
	  if (DBL_LT(0, xScaledTemp(i))) {
	    xScaled(i) = xScaledTemp(i);
	  }
	}
	
	manager.computeSmoothComponent(false /* no need for hessian */);
	double f2  = manager.getSmoothObj();
	const AccpmVector *df2 = manager.getSmoothGradient(true);
	//r = (xScaled' * LocSetS.AT)' + df2;
	AccpmVector r(my);
	AccpmLAMatVecMult(locSet.getA(), xScaled, r, 1, 0);
	if (df2) {
	  AccpmLAAddMult(r, 1, *df2);
	}

	AccpmVector dinf;
	AccpmVector dsup;

	if (_param->getIntParameter("Box")) {
	  const AccpmVector *varB = manager.getVariableLB(true);
	  if (varB) {
	    dinf = y;
	    AccpmLAAddMult(dinf, -1, *varB);
	    dinf.invert();
	    AccpmLAScale(scaleFactor, dinf);
	    AccpmLAAddMult(r, -1, dinf);
	  }
	  varB = manager.getVariableUB(true);
	  if (varB) {
	    dsup = *varB;
	    AccpmLAAddMult(dsup, -1, y);
	    dsup.invert();
	    AccpmLAScale(scaleFactor, dsup);
	    AccpmLAAddMult(r, 1, dsup);
	  }
	}

	lBound = f2 - AccpmLADotProd(locSet.getC(), xScaled) + AccpmLADotProd(r, y);
	if (df2) {
	  lBound -= AccpmLADotProd(*df2, y);
	}
	AccpmVector centeredY = y;
	AccpmLAAddMult(centeredY, -1, manager.getBestY(true));
	lBound -= _param->getRealParameter("Delta") * AccpmLANorm2(centeredY) * AccpmLANorm2(r);
	
	if (_param->getIntParameter("Box")) {
	  const AccpmVector *varB = manager.getVariableLB(true);
	  if (varB) {
	    lBound += AccpmLADotProd(*varB, dinf);
	  }
	  varB = manager.getVariableUB(true);
	  if (varB) {
	    lBound -= AccpmLADotProd(*varB, dsup);
	  }
	}
	
      }
    }
    manager.getPointGen()._xScaled = xScaled;
  }
  if (_computeLB) {
    if (IS_FINITE(lBound) && (DBL_LT(manager.getObjUB(), lBound))) {
      AccpmWarning("Lower bound was overestimated : use previous value");
      lBound = ACCPM_MINUS_INF;
    }
    manager.updateLB(lBound);
  }
  return lBound;
}

/**
 * This function solves the Schur complement system $(a - b * b') x = rhs$
 * via Sherman-Morrison formula. The matrix [ a b ; b' 1 ] is assumed
 * to be positive definite.
 */
void
DualMethod::solveSCSystem(const AccpmVector &a, const AccpmVector &b, const AccpmVector &rhs,
                          AccpmVector &x) const
{
  AccpmVector bBya = b;
  bBya.rdivide(a);
  double sc = 1 - AccpmLADotProd(b, bBya);
  AccpmVector x1 = rhs;
  x1.rdivide(a);
  double scale = AccpmLADotProd(b, x1)/sc;
  x = x1;
  AccpmLAAddMult(x, scale, bBya);
}

void
DualMethod::handleEqualityConstraints(Manager &manager, LocSet &locSet, 
				      const AccpmVector &wsV, double ws0, 
				      double f2, const AccpmVector *df2,
				      const AccpmGenMatrix *d2f2,
				      const AccpmVector &diagQ,
				      const NewtonSolution &sol,
				      NewtonSolution &step) const
{
  AccpmGenMatrix *constraints = 0;
  AccpmVector *rhs = 0;
  _param->getEqualityConstraints(constraints, rhs);
  if (constraints) {
    assert(rhs);
    int my = step._y.size();
    int n = step._s.size();
    int p = step._z.size();
    int numConstraints = rhs->size();
    AccpmVector rd(numConstraints);
    AccpmLAMatTransVecMult(*constraints, sol._y, rd, -1, 0);
    AccpmLAAddMult(rd, 1, *rhs);
    AccpmGenMatrix Ny(my, my);
    AccpmGenMatrix Ns(n, my);
    AccpmGenMatrix Nz;
    AccpmVector Ns0;
    AccpmVector Nss;
    AccpmVector Nzs;
    if (manager.inPhase2()) {
      Nz.resize(p, my);
    }
    AccpmVector ry(my);
    ry = 0;
    AccpmVector rz(p);
    rz = 0;
    AccpmVector rs(n);
    rs = 0;
    for (int i = 0; i < my; ++i) {
      ry(i) = -1;
      NewtonSolution tmpstep(my, n, p);
      newtonDualDir(manager, locSet, ry, rz, rs, 0, wsV, ws0, f2, df2, d2f2, diagQ, sol, tmpstep, true);
      Ny.assignColumn(i, tmpstep._y);
      Ns.assignColumn(i, tmpstep._s);
      if (manager.inPhase2()) {
	Nz.assignColumn(i, tmpstep._z);
	Ns0.append(tmpstep._s0);
	Nss.append(tmpstep._ss);
	Nzs.append(tmpstep._zs);
      }
      ry(i) = 0;
    }
    AccpmVector eX(numConstraints);
    AccpmGenMatrix DND(numConstraints, numConstraints);
    AccpmGenMatrix tmp(my, numConstraints);
    AccpmLAMatMatMult(Ny, *constraints, tmp, 1, 0);
    AccpmLAMatTransMatMult(*constraints, tmp, DND, 1, 0);
    
    // X = (DtnyD)^-1(rd - Dtdy1)
    AccpmLAMatTransVecMult(*constraints, step._y, rd, -1, 1);
    AccpmLALinearSolve(DND, eX, rd);
    AccpmVector Dx(my);
    AccpmLAMatVecMult(*constraints, eX, Dx, 1, 0);
    NewtonSolution step2(my, n, p);
    AccpmLAMatVecMult(Ny, Dx, step2._y, 1, 0);
    AccpmLAMatVecMult(Ns, Dx, step2._s, 1, 0);
    if (manager.inPhase2()) {
      AccpmLAMatVecMult(Nz, Dx, step2._z, 1, 0);
      step2._s0 = AccpmLADotProd(Ns0, Dx);
      step2._ss = AccpmLADotProd(Nss, Dx);
      step2._zs = AccpmLADotProd(Nzs, Dx);
    } else {
      step2._z = 0;
    }
    step.addMult(step2, 1);
  }

}

int DualMethod::linSolve(Manager &manager, const LocSet &locSet,
			 const AccpmGenMatrix &A, RealMatrix &X,
			 const RealMatrix &B) const
{
  int status = AccpmLASymmLinSolve(A, X, B);
  if (status != 0) {
    manager.setExitCode(LA_ERROR);
  }
  AccpmVector ky;
  AccpmVector kx;
  if (status > 0) {
    if (_param->getIntParameter("CheckLocSetInterior")) {
      double kellyBound;
      int feasible = locSet.checkFeasibility(*_param, kellyBound, ky, kx);
      if (feasible != 1) {
	status = 1;
      } else if (feasible == 1) {
	std::cout << "Updating lower bound with the obtained Kelly bound: " 
		  << kellyBound << std::endl;
	manager.updateLB(kellyBound);
	manager.getPointGen()._xScaled = kx;
	status = 1;
      }
    }
  }
  return status;

}

}
