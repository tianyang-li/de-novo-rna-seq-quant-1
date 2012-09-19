// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef SOLUTION_H
#define SOLUTION_H

/**
 * @file Method.h
 * @ingroup AccpmCore
 */
#include "AccpmVector.h"

namespace Accpm 
{
  
/**
 * Abstract class for solution for the methods
 *
 * @ingroup AccpmCore */  
  class Solution {
  public:
    Solution() {};
    virtual ~Solution() {};
  };
  class Manager;  
/**
 * Class for keepin the solution of the Newton system
 *
 * @ingroup AccpmCore */  
 
  class NewtonSolution : public Solution {
  
  public:
    NewtonSolution(int m = 0, int n = 0, int p = 0);
    ~NewtonSolution();

    void init(const Manager &manager);
    int warmStart(const Manager &manager);
    void output(std::ostream &os) const;
    void addMult(const NewtonSolution &sol, double alpha = 1.0);
						
    AccpmVector _y;
    AccpmVector _s;
    AccpmVector _z;
    double _s0;
    double _zs;
    double _ss;
    
  };
}
#endif
