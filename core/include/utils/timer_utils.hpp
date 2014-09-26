/****************************************************************************************
 * File: timer_utils.hpp
 * Description: This module provides a tool for getting elapsed time.
 * 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 30/12/2013
 ***************************************************************************************/

#ifndef UTILS_TIMER_UTILS_HPP
#define UTILS_TIMER_UTILS_HPP


#ifdef WIN32   // Windows system specific
#include <windows.h>
#else      // Unix based system specific
#include <sys/time.h>
#endif
#include <sstream>

#include <stdlib.h>
#include <string>

namespace utility 
{

class Timer
{
 public:
  Timer() {
#ifdef WIN32
    QueryPerformanceFrequency(&frequency);
    startCount.QuadPart = 0;
    endCount.QuadPart = 0;
#else
    startCount.tv_sec = startCount.tv_usec = 0;
    endCount.tv_sec = endCount.tv_usec = 0;
#endif
    stopped = 0;
    startTimeInMicroSec = 0;
    endTimeInMicroSec = 0;
  }
  ~Timer() {}                                   // default destructor
  
   void start();                             // starts the timer.
   void stop();                              // stops the timer.
   double getElapsedTime();                    // gets elapsed time in seconds.
   double getElapsedTimeInSec();               // gets elapsed time in seconds (same as getElapsedTime).
   double getElapsedTimeInMilliSec();          // gets elapsed time in milli-seconds.
   double getElapsedTimeInMicroSec();          // gets elapsed time in micro-seconds.
   void restart();

  std::string display();
  
 protected:
 private:
  double startTimeInMicroSec;                 // starting time in micro-seconds
  double endTimeInMicroSec;                   // ending time in micro-seconds
  int    stopped;                             // stop flag 
#ifdef WIN32
  LARGE_INTEGER frequency;                    // ticks per second
  LARGE_INTEGER startCount;                   
  LARGE_INTEGER endCount;
#else
  timeval startCount;                         
  timeval endCount;                           
#endif
};


inline std::string timeDisplay( double seconds ) {
  std::ostringstream result;
  long secs = seconds / 1L; // 1L: a long with value 1
  int hours = (secs / 3600);
  int minutes = (secs / 60) % 60;
  int s = secs % 60;

  if ( secs < 60 )  {
    result << seconds << " (seconds).";

  } else if (secs < 3600) {
    result << minutes << " (minutes) " << s << " (seconds).";

  } else {
    result << hours << " (hours), " << minutes << " (minutes) " << s << " (seconds)";
  }  

  return result.str();
}


} //namespace utils ends here.


#endif // UTILS_TIMER_UTILS_HPP
