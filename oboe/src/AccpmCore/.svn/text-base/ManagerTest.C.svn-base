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

#include "Manager.h"
#include "Parameters.h"
#include "Method.h"
#include "LocSet.h"
#include <iostream>

int 
main(int argc, char *argv[])
{
  // Need to create Parameter before
  Accpm::Parameters param("param.txt");
  std::cout << "Parameters:\n" <<  param << std::endl;
  std::cout << "Creating Manager" << std::endl;
  Accpm::Manager manager(&param);
  manager.check();
  Accpm::DualMethod dMethod(&param);
//  dMethod.run(manager);
  Accpm::LocSet locSet(manager, param);

}
