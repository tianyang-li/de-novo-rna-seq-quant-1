// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

/**
 * @file Method.h
 * @ingroup AccpmCore
 */
#include "Solution.h"
#include "Manager.h"
#include "AccpmDefs.h"

#include "AccpmBlasInterface.h"

namespace Accpm {

NewtonSolution::NewtonSolution(int m, int n, int p)  : Solution(), _y(m), _s(n), _z(p), _s0(0), _zs(0), _ss(0)
{
  _y = 0;
  _s = 0;
  _z = 0;
}

NewtonSolution::~NewtonSolution() 
{
}

void
NewtonSolution::init(const Manager &manager) 
{
  _y = manager.getActiveY();
  _z = manager.getCurrentZ();
  warmStart(manager);

  if (manager.inPhase2()) {
    _ss = manager.getPointGen()._ss;
    _zs = manager.getPointGen()._zs;
    /* double f2 = manager.getSmoothObj();
    _zs = f2;//std::max(f2*(1 + 1e-10), 1.0);
    _ss = 100;//_zs - f2;
    _s0 = 100;  */
  } else {
    _ss = manager.getPointGen()._ss;
    _zs = manager.getPointGen()._zs;
  }
}

int
NewtonSolution::warmStart(const Manager &manager)
{ 
  if (manager.getSizeOfs() == 0) {
    _s = 1;
    if (manager.inPhase2()) {
      _s0 = 1;
    }
    return 0;
  } 
  int n = _s.size();
  AccpmVector sp;
  double smin = ACCPM_PLUS_INF;
  int ls = manager.getPointGen()._s.size();
  if (ls > 1) {
  
    for (int i = 0; i < ls - 1; ++i) {
      if (smin > manager.getPointGen()._s(i)) {
	smin = manager.getPointGen()._s(i);
      }
    }
    if (manager.inPhase2() && manager.previousPhase2()) {
      sp = AccpmVector(n - ls + 1);
      sp = smin;
    } else {
      sp = AccpmVector(n - ls);
      if (smin > manager.getPointGen()._s(ls-1)) {
	smin = manager.getPointGen()._s(ls-1);
      }
      sp = smin;
    }
  } else {
    sp = AccpmVector(n-ls);
    sp = 1;
  }
  if (manager.inPhase2()) {
    if (manager.previousPhase2()) {
      assert(n == ls-1 + sp.size());
      assert(ls >= 2);
      _s = AccpmVector(n);
      LaIndex index(0, ls-2);
	  // Following inject causes trouble on WINDOWS .NET 2003
     // _s(index).inject(manager.getPointGen()._s(index));
      for (int i = 0; i < ls-1; ++i) {
        _s(i) = manager.getPointGen()._s(i);
	  }
	  for (int i = ls - 1; i < n; ++i) {
	_s(i) = sp(i-ls+1);
      }
      _s0 = manager.getPointGen()._s(ls-1);
    } else {
      _s = manager.getPointGen()._s;
      _s.append(sp);
      _s0 = sp(0);
    }
  } else {
    if (ls >= 1) {
      _s = manager.getPointGen()._s;
      _s.append(sp);
    } else {
      _s = sp;
    }
  }
  return 0;
}


void 
NewtonSolution::output(std::ostream &os) const
{
  os << "y: \n" << _y << std::endl;
  os << "z: \n" << _z << std::endl;
  os << "s: \n" << _s << std::endl;
  os << "s0: " << _s0 << std::endl;
  os << "ss: " << _ss << std::endl;
  os << "zs: " << _zs << std::endl;
}

void
NewtonSolution::addMult(const NewtonSolution &sol, double alpha)
{
  AccpmLAAddMult(_y, alpha, sol._y);
  AccpmLAAddMult(_s, alpha, sol._s);
  AccpmLAAddMult(_z, alpha, sol._z);
  _s0 += alpha * sol._s0;
  _ss += alpha * sol._ss;
  _zs += alpha * sol._zs;
}

}
