// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "LocSet.h"
#include "Parameters.h"

#include <cmath>
#include "AccpmBlasInterface.h"

namespace Accpm {

LocSet::LocSet(const Manager &manager, const Parameters &param) 
{
  manager.getActiveCuts(_A, true);
  _proximalCenter = manager.getBestY(true);
  
  _c = manager.getRhsCoef();
  if (manager.inPhase2()) {
    int numCuts = _A->size(1);
    int p = param.getIntParameter("NumSubProblems");
    _E = AccpmGenMatrix(p, numCuts);
    _E = 0;
    const AccpmVector &subProblemIndex = manager.getSubProblemIndex();
    for (int i = 0; i < numCuts; ++i) {
      int index = (int)subProblemIndex(i);
      if (index > 0) {
	_E(index-1, i) = 1;
      }
    }
    //std::cout << "E:\n" << _E << std::endl;
    //std::cout << "c:\n" << _c << std::endl;
    _rhs = manager.getObjUB(); 
    if (IS_FINITE(_rhs)) {
      _rhs += 1e-5 * (param.getRealParameter("Tolerance")) * std::fabs(_rhs);
    }
  }
  _AT = _A->transpose();
  _AFull = 0;
  _EFull = 0;
  _ATQAFull = 0;
  _EFullT = 0;
}

LocSet::~LocSet()
{
  clear();
 
  if (_AT) {
    delete _AT;
    _AT = 0;
  }
  if (_ATQAFull) {
    delete _ATQAFull;
    _ATQAFull = 0;
  }
}

void
LocSet::clear()
{ 
  delete _AFull;
  _AFull = 0;
  
  delete _EFull;
  _EFull = 0;
  
  delete _EFullT;
    _EFullT = 0;
}

int 
LocSet::computeFullAE(const Parameters &param, const AccpmVector *df2)
{

  // LocSetS.AFull = [LocSetS.A zeros(my,1) df2];
  clear();
  int my = _A->size(0);
  int n = _A->size(1);
  int p = _E.size(0);
  _AFull = new AccpmGenMatrix(my, n + 2);
  (*_AFull)(_A->index(0), _A->index(1)).inject(*_A);
  AccpmVector v(my);
  v = 0;
  _AFull->assignColumn(n, v);
  if (df2) {
    _AFull->assignColumn(n+1, *df2);
  } else {
    _AFull->assignColumn(n+1, v);
  }
  // LocSetS.EFull = [LocSetS.E  -LocSetS.pi zeros(mz,1); zeros(1,n) -1 1];
  _EFull = new AccpmGenMatrix(p+1, n+2);
  *_EFull = 0;
  (*_EFull)(_E.index(0), _E.index(1)).inject(_E);
  AccpmVector pi;
  pi = *param.getPi();
  pi.negate();
  pi.append(-1);
  _EFull->assignColumn(n, pi);
  (*_EFull)(p, n+1) = 1;
  //std::cout << "AFull\n" << *_AFull << std::endl;
  //std::cout << "EFull\n" << *_EFull << std::endl;
  
  _EFullT = _EFull->transpose();
  return 0;
}

int 
LocSet::computeFullATQA(const AccpmGenMatrix &ATQA, const AccpmVector &diagQ, const AccpmVector *df2)
{
  int my = _A->size(0);
  int n = _A->size(1);
  if (_ATQAFull == 0) {
    _ATQAFull = new AccpmGenMatrix(n+2, n+2);
    *_ATQAFull = 0;
  }

  if (df2) {
    // bTQ =  [zeros(my,1)./ diaQ  df2 ./ diaQ];
    AccpmGenMatrix bTQ(my, 2);
    bTQ = 0;
    // bTQb = bTQ' * [zeros(my,1) df2];
    AccpmGenMatrix bTQb(2,2);
    bTQb = 0;
    // bTQA = bTQ' * LocSetS.A;
    AccpmGenMatrix bTQA(2, n);
    bTQA = 0;
    AccpmGenMatrix bTQAT(n, 2);
    bTQA = 0;
    AccpmVector v = *df2;
    v.rdivide(diagQ);
    bTQ.assignColumn(1, v);
    bTQb(1,1) = AccpmLADotProd(v, *df2);
    
    AccpmVector result(n);
    AccpmLAMatTransVecMult(*_A, v, result, 1, 0);
    bTQA.assignRow(1, result);
    bTQAT.assignColumn(1, result);
    
    LaIndex rowIndex(n, n+1);
    (*_ATQAFull)(rowIndex, ATQA.index(1)).inject(bTQA);
    (*_ATQAFull)(ATQA.index(0), rowIndex).inject(bTQAT);
    (*_ATQAFull)(rowIndex, rowIndex).inject(bTQb);
  }
  
  // LocSetS.ATQAFull = Concat(PointGenS.ATQA, bTQA', bTQA, bTQb);
  (*_ATQAFull)(ATQA.index(0), ATQA.index(1)).inject(ATQA);
  //std::cout << "ATQAFull\n" << *_ATQAFull << std::endl;

  return 0;
}

int
LocSet::checkFeasibility(const Parameters &param, double &objVal, 
			 AccpmVector &y, AccpmVector &x) const
{
  if (!param.getIntParameter("CheckLocSetInterior")) {
    return 1;
  }
#ifdef OBOE_HAS_GLPK
  return solveLP(param, objVal, y, x);
#else
#ifdef COIN_HAS_OSI
  return solveLP(param, objVal, y, x);
#else
  AccpmError("No LP Solver present. Configure --with-glpk and try again.");
  objVal = ACCPM_MINUS_INF;
  return 1;
#endif
#endif

}

int
LocSet::solveLP(const Parameters &param, double &objVal, 
		      AccpmVector &y, AccpmVector &x) const
{
  
  const string lpSolver = param.getStringParameter("LPSolverName");
  int status = 0;

#ifdef COIN_HAS_OSI
  OsiSolverInterface *si;
  if (lpSolver == "GLPK") {
    si = new OsiGlpkSolverInterface;
  } else if (lpSolver == "CLP") {
    si = new OsiClpSolverInterface;
  } else {
    AccpmError("Only GLPK or CLP LP solver supported.");
    exit(0);
  }

  status = solveWithOsi(si, param, objVal, y, x);
  delete si;
  return status;
#else
#ifdef OBOE_HAS_GLPK
  status = solveWithGLPK(param, objVal, y, x);
#endif
#endif
  return status;
}
#ifdef COIN_HAS_OSI
int 
LocSet::solveWithOsi(OsiSolverInterface *si, const Parameters &param, 
		     double &objVal, AccpmVector &y, AccpmVector &x) const
{
#ifdef COIN_HAS_OSI
 
 int nrows = _A->size(1);
 assert(nrows == _c.size());
 int offset = _A->size(0) + _E.size(0);
 int ncols = offset + nrows;
 double plusInf = si->getInfinity();
 double minusInf = -plusInf;
 double *row_lb = new double[nrows]; //the row lower bounds
 double *row_ub = new double[nrows]; //the row upper bounds        
 for (int i = 0; i < nrows; ++i) {
   //lpx_set_row_bnds(lp, i+1, LPX_UP, 0, _c(i));
   row_ub[i] = _c(i);
   row_lb[i] = minusInf;
 }
 
 const AccpmVector *b = param.getB();
 const AccpmVector *lb = param.getVariableLB();
 const AccpmVector *ub = param.getVariableUB();
 
 double *objective    = new double[ncols];//the objective coefficients
 memset(objective, 0, sizeof(double) * ncols);
 
 double *col_lb       = new double[ncols];//the column lower bounds
 memset(col_lb, minusInf, sizeof(double) * ncols);
 
 double *col_ub       = new double[ncols];//the column upper bounds
 memset(col_ub, plusInf, sizeof(double) * ncols);

 if (lb && ub) {
    for (int i = 0; i < _A->size(0); ++i) {
      col_lb[i] = (*lb)(i);
      col_ub[i] = (*ub)(i);
      //  lpx_set_col_bnds(lp, i+1, LPX_DB, (*lb)(i), (*ub)(i));
    }
  } else if (lb) {
    for (int i = 0; i < _A->size(0); ++i) {
      col_lb[i] = (*lb)(i);
      col_ub[i] = plusInf;
      // lpx_set_col_bnds(lp, i+1, LPX_LO, (*lb)(i), 0);
    }
  } else if (ub) {
    for (int i = 0; i < _A->size(0); ++i) {
      col_lb[i] = minusInf;
      col_ub[i] = (*ub)(i); 
      //lpx_set_col_bnds(lp, i+1, LPX_UP, 0, (*ub)(i));
    }
  } else { // No bounds
    for (int i = 0; i < _A->size(0); ++i) {
      col_lb[i] = minusInf;
      col_ub[i] = plusInf;
      //lpx_set_col_bnds(lp, i+1, LPX_FR, 0, 0);
    }
  }
  if (b) {
    if (param.getOptType() == OPT_MAX) {
      for (int i = 0; i < _A->size(0); ++i) {
	objective[i] = -(*b)(i);
	//lpx_set_obj_coef(lp, i+1, -(*b)(i));
      }
    } else {
      for (int i = 0; i < _A->size(0); ++i) {
	objective[i] = (*b)(i);
	//lpx_set_obj_coef(lp, i+1, (*b)(i));
      }
    }
  }

 /* The z variables are free */
  AccpmVector pi(_E.size(0));
  if (param.getPi()) {
    pi = *param.getPi();
  } else {
    pi = 1;
  }
  for (int i = _A->size(0), j = 0; i < offset; ++i, ++j) {
    col_lb[i] = minusInf;
    col_ub[i] = plusInf;
    //    lpx_set_col_bnds(lp, i+1, LPX_FR, 0, 0);
    objective[i] = pi(j);
    //    lpx_set_obj_coef(lp, i+1, pi(j));
  }
    
  for (int i = 0; i < nrows; ++i) {
    /* The slacks are bounded from below by 0 */
    col_lb[offset+i] = 0;
    col_ub[offset+i] = plusInf;
    //lpx_set_col_bnds(lp, offset+i+1, LPX_LO, 0, 0);

    /* The slacks have an objective coefficient, which is 1 */
    objective[offset+i] = 1;
    //lpx_set_obj_coef(lp, offset+i+1, 1);
  }
 
  /* Add A part of the constraint matrix Aty*/
  int *ia = new int [nrows*ncols+1];
  int *ja = new int [nrows*ncols+1];
  double *ar = new double [nrows*ncols+1];
  for (int i = 0; i < nrows*ncols+1; ++i) {
    ia[i] = 0;
    ja[i] = 0;
    ar[i] = 0;
  }

  int nz = 0;
  for (int j = 0; j < _AT->size(1); ++j) {
    AccpmVector col = _AT->getColumn(j);
    for (int i = 0; i < col.size(); ++i) {
      if (!DBL_CMP(col(i), 0)) {
	++nz;
	ia[nz] = i;
	ja[nz] = j;
	ar[nz] = col(i);
      }
    }
  }
offset = _A->size(0);
  /* Add E part of the constraint matrix -Etz*/  
  for (int i = 0; i < _E.size(1); ++i) {
    AccpmVector row = _E.getColumn(i);
    for (int j = 0; j < row.size(); ++j) {
      if (!DBL_CMP(row(j), 0)) {
	++nz;
	ia[nz] = i;
	ja[nz] = offset+j;
	ar[nz] = -row(j);
      }
    }
  }
  
  offset = _A->size(0) + _E.size(0);
  /* Add the slack columns */
  for (int i = 0; i < nrows; ++i) {
    ++nz;
    ia[nz] = i;
    ja[nz] = offset+i;
    ar[nz] = -1; /* Aty - Etz - s <= c, so positive slacks => infeasible */
  }
 
  CoinPackedMatrix matrix(true, ia, ja, ar, nz);
  si->loadProblem(matrix, col_lb, col_ub, objective, row_lb, row_ub);
  si->writeMps("locset");
  si->initialSolve();
  int feasible = 1;
  
  if ( si->isProvenOptimal() ) {
    std::cout << "Found optimal solution!" << std::endl;
    objVal =  si->getObjValue();
    std::cout << "Objective value is " << objVal << std::endl;
    
    //int n = si->getNumCols();
    const double *solution = si->getColSolution();
    const double *dualVar = si->getRowPrice();
    x.resize(nrows, 1);
    for (int i = 0; i < nrows; ++i) {
      //double s = lpx_get_col_prim(lp, offset+i+1);
      double s = solution[offset+i];
      //x(i) = -lpx_get_row_dual(lp, i+1);
      x(i) = -dualVar[i];
      if (DBL_GT(s, 0)) {
	feasible = 0;
	if (param.getIntParameter("Verbosity") > 1) {
	  std::cout << "Constraint " << i << " is not feasible: infeasibility: " 
		    << s << std::endl;
	}
      }
    }
    y.resize(_A->size(0), 1);
    for (int i = 0; i < _A->size(0); ++i) {
      y(i) = solution[i];
	// y = lpx_get_col_prim(lp, i+1);
    }
    //  std::cout << "x:" << x << "\nx.sum()" << x.sum() 
    //	    << " cx: " << AccpmLADotProd(x, _c) << std::endl;
    if (!feasible) {
      objVal = ACCPM_MINUS_INF;
      std::cout << "Localization set is infeasible" << std::endl;
    } else {
      std::cout << "Localization set is feasible" << std::endl;
    }
    
    delete ia;
    delete ja;
    delete ar;

    delete row_lb;
    delete row_ub;
    delete col_lb;
    delete col_ub;
    delete objective;

  return feasible;
    // We could then print the solution or examine it.
  } else {
    std::cout << "Didn't find optimal solution." << std::endl;
    // Could then check other status functions.
    objVal = ACCPM_MINUS_INF;
    return 0;
  }
   
#else 
  return 1;
#endif

}
#endif

int 
LocSet::solveWithGLPK(const Parameters &param, double &objVal, 
		      AccpmVector &y, AccpmVector &x) const
{
#ifdef OBOE_HAS_GLPK
  LPX *lp = lpx_create_prob();
  lpx_set_prob_name(lp, "localization set");
  lpx_set_obj_dir(lp, LPX_MIN);
  int nrows = _A->size(1);
  lpx_add_rows(lp, nrows);
  assert(nrows == _c.size());
  for (int i = 0; i < nrows; ++i) {
    lpx_set_row_bnds(lp, i+1, LPX_UP, 0, _c(i));
  }

  int offset = _A->size(0) + _E.size(0);
  int ncols = offset + nrows;
  const AccpmVector *b = param.getB();
  const AccpmVector *lb = param.getVariableLB();
  const AccpmVector *ub = param.getVariableUB();

  lpx_add_cols(lp, ncols);

  /* Set the bounds on the y variables */
  if (lb && ub) {
    for (int i = 0; i < _A->size(0); ++i) {
      lpx_set_col_bnds(lp, i+1, LPX_DB, (*lb)(i), (*ub)(i));
    }
  } else if (lb) {
    for (int i = 0; i < _A->size(0); ++i) {
      lpx_set_col_bnds(lp, i+1, LPX_LO, (*lb)(i), 0);
    }
  } else if (ub) {
    for (int i = 0; i < _A->size(0); ++i) {
      lpx_set_col_bnds(lp, i+1, LPX_UP, 0, (*ub)(i));
    }
  } else { // No bounds
    for (int i = 0; i < _A->size(0); ++i) {
      lpx_set_col_bnds(lp, i+1, LPX_FR, 0, 0);
    }
  }
  if (b) {
    if (param.getOptType() == OPT_MAX) {
      for (int i = 0; i < _A->size(0); ++i) {
	lpx_set_obj_coef(lp, i+1, -(*b)(i));
      }
    } else {
      for (int i = 0; i < _A->size(0); ++i) {
	lpx_set_obj_coef(lp, i+1, (*b)(i));
      }
    }
  }
  
  /* The z variables are free */
  AccpmVector pi(_E.size(0));
  if (param.getPi()) {
    pi = *param.getPi();
  } else {
    pi = 1;
  }
  for (int i = _A->size(0), j = 0; i < offset; ++i, ++j) {
    lpx_set_col_bnds(lp, i+1, LPX_FR, 0, 0);
    lpx_set_obj_coef(lp, i+1, pi(j));
  }
    
  for (int i = 0; i < nrows; ++i) {
    /* The slacks are bounded from below by 0 */
    lpx_set_col_bnds(lp, offset+i+1, LPX_LO, 0, 0);
    /* The slacks have an objective coefficient, which is 1 */
    lpx_set_obj_coef(lp, offset+i+1, 1);
  }

 
  /* Add A part of the constraint matrix Aty*/
  int *ia = new int [nrows*ncols+1];
  int *ja = new int [nrows*ncols+1];
  double *ar = new double [nrows*ncols+1];
  for (int i = 0; i < nrows*ncols+1; ++i) {
    ia[i] = 0;
    ja[i] = 0;
    ar[i] = 0;
  }

  int nz = 0;
  for (int j = 0; j < _AT->size(1); ++j) {
    AccpmVector col = _AT->getColumn(j);
    for (int i = 0; i < col.size(); ++i) {
      if (!DBL_CMP(col(i), 0)) {
	++nz;
	ia[nz] = i+1;
	ja[nz] = j+1;
	ar[nz] = col(i);
      }
    }
  }

  /*
    std::cout << "A nzeroes: " << nz
    <<" nrows: " << nrows << " ncols: " << _A->size(0) << std::endl;
  */

  offset = _A->size(0);
  /* Add E part of the constraint matrix -Etz*/  
  for (int i = 0; i < _E.size(1); ++i) {
    AccpmVector row = _E.getColumn(i);
    for (int j = 0; j < row.size(); ++j) {
      if (!DBL_CMP(row(j), 0)) {
	++nz;
	ia[nz] = i+1;
	ja[nz] = offset+j+1;
	ar[nz] = -row(j);
      }
    }
  }
  
  offset = _A->size(0) + _E.size(0);
  /* Add the slack columns */
  for (int i = 0; i < nrows; ++i) {
    ++nz;
    ia[nz] = i+1;
    ja[nz] = offset+i+1;
    ar[nz] = -1; /* Aty - Etz - s <= c, so positive slacks => infeasible */
  }
 

  lpx_load_matrix(lp, nz, ia, ja, ar);

  lpx_set_int_parm(lp, LPX_K_PRESOL, 1);
  /* Solve the LP */
  if (lpx_simplex(lp) != LPX_E_OK) {
    return -1;
  }
  
  int feasible = 1;
  x.resize(nrows, 1);
  for (int i = 0; i < nrows; ++i) {
    double s = lpx_get_col_prim(lp, offset+i+1);
    x(i) = -lpx_get_row_dual(lp, i+1);
    if (DBL_GT(s, 0)) {
      feasible = 0;
      if (param.getIntParameter("Verbosity") > 1) {
	std::cout << "Constraint " << i << " is not feasible: infeasibility: " 
		  << s << std::endl;
      }
    }
  }
  y.resize(_A->size(0), 1);
  for (int i = 0; i < _A->size(0); ++i) {
    y(i) = lpx_get_col_prim(lp, i+1);
  }
  //  std::cout << "x:" << x << "\nx.sum()" << x.sum() 
  //	    << " cx: " << AccpmLADotProd(x, _c) << std::endl;
  if (!feasible) {
    objVal = ACCPM_MINUS_INF;
    std::cout << "Localization set is infeasible" << std::endl;
  } else {
    objVal = lpx_get_obj_val(lp);
    std::cout << "Localization set is feasible" << std::endl;
  }

  delete ia;
  delete ja;
  delete ar;
  lpx_delete_prob(lp);

  return feasible;

#else 
  return 1;
#endif

}

}
