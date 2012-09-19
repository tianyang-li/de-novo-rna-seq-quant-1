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

#include "Parameters.h"
#include <iostream>

int 
main(int argc, char *argv[])
{

  Accpm::Parameters param("param.oboe");
  std::cout << "Parameters Read from file" << param << std::endl;
  std::cout << "Getting Num Var: ";
  int i = param.getIntParameter("NumVariables");
  std::cout << i << std::endl;
  param.setIntParameter("NumVariables", 10);
  i = param.getIntParameter("NumVariables");
  std::cout << "Modified value: " << i << std::endl;
  assert(i == 10);

  param.setStringParameter("MethodName", "Foo");
  
  vector<double> start(10, 1);
  param.setStartingPoint(start);
  param.setPi(start);
  param.setIntParameter("Verbosity", 2);
  param.setOptimizationType("Max");
  param.setB(start);
  param.setIntParameter("Ball", 1);
  param.setRealParameter("RadiusBall", 10);
  std::cout << "Modified Parameters" << param << std::endl;

 
}
