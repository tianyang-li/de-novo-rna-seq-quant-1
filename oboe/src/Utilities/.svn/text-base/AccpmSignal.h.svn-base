// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// The OBOE team
//

#ifndef ACCPM_SIGNAL_H
#define ACCPM_SIGNAL_H

#include <fenv.h>
#include <signal.h>

static void 
signalHandler(int signum)
{
  std::cout << "Signal Handler called with signal " 
	    << signum << std::endl;
  //longjmp(context, signum);
  exit(0);
}

void
AccpmInstallSignalHandler()
{
  feenableexcept(FE_UNDERFLOW);
  signal(SIGFPE, signalHandler);
}

#endif
