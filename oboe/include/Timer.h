// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// The OBOE team
//

#ifndef TIMER_H
#define TIMER_H

#include <config.h>
#include <iostream>
#if (defined(LINUX))
#include <sys/times.h>
#include <unistd.h>
#endif

#if (defined(WIN32))
#include <time.h>
//#include <unistd.h>
#endif

/**
 * The Timer class. 
 * Takes care of tracking Real, CPU and System time.
 **/  

namespace Accpm {

class Timer {

#if (defined(LINUX))
 private:
  struct tms _startTime;
  clock_t _startClock;
  struct tms _endTime;
  clock_t _endClock;
  const double _clockTicks;

 public:
  Timer() : _clockTicks(sysconf(_SC_CLK_TCK)) {}
    inline void start() { _startClock = times(&_startTime); }
    inline void stop()  { _endClock = times(&_endTime); }
    double getRealTime() const { return (_endClock - _startClock)/_clockTicks;}
    double getCpuTime()  const { return (_endTime.tms_utime - _startTime.tms_utime 
					 + _endTime.tms_stime - _startTime.tms_stime)/_clockTicks; }
    double getSysTime()  const { return (_endTime.tms_stime - _startTime.tms_stime)/_clockTicks; }
#endif

#if (defined(WIN32))
 private:
  time_t _startTime;
  clock_t _startClock;
  time_t _endTime;
  clock_t _endClock;
  const double _clockTicks;
  
 public:
  Timer() : _clockTicks(CLK_TCK) {}
    inline void start() { _startClock = time(&_startTime); _startClock = clock_t();}
    inline void stop()  { _endClock = time(&_endTime); _endClock = clock_t();}
    double getRealTime() const { return _endTime - _startTime;}
    //double getCpuTime()  const { return  (_endClock - _startClock)/_clockTicks; }
	// Dont know how to get CPU time on windows
	double getCpuTime()  const { return  getRealTime(); }
    double getSysTime()  const { std::cout << "Function Timer::getSysTime() not supported\n"; return 0; }
#endif

};


}
#endif
