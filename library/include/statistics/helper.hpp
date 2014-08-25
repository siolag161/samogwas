

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

/*** ***/

#ifndef __HELPER_H__
#define __HELPER_H__

#include <string>
#include <vector>
#include <cstdio>
#include <iostream>
#include <sstream>

//#include "plink.h"
//#include "options.h"

namespace stats {

template<class T>
inline const T SQR(const T a) {return a*a;}
 
template<class T>
inline const T MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}
 
template<class T>
inline const T MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}
 
template<class T>
inline const T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
 
template<class T>
inline void SWAP(T &a, T &b)
{T dum=a; a=b; b=dum;}


inline void error(std::string msg);

inline std::string int2str(int);
inline std::string dbl2str(double,int prc = -1);
inline std::string dbl2str_fixed(double, int prc = -1);
inline std::string longint2str(long int);
inline std::string sw(std::string s , int n);
inline std::string sw(double d , int n);
inline std::string sw(double d , int f, int n);
inline std::string sw(int i , int n);

inline std::string itoa(int, int);

inline bool realnum(double d);
}
/////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <ctime>

namespace stats {

inline void error(std::string msg)
{
  std::cerr << "\nERROR: " << msg << "\n";
  // LOG << "\nERROR: " << msg << "\n";  
  // LOG.close();

  // if (par::gplink)
  //   {
  //     ofstream GP((par::output_file_name+".gplink").c_str(),ios::out);
  //     GP << "1\n";
  //     GP.close();
  //   }
  
  // PP->cleanUp();

  exit(1);
}

inline std::string sw(std::string s , int n) 
{
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

inline std::string sw(double d , int n) 
{
  std::string s = realnum(d) ? dbl2str(d) : "NA";
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

inline std::string sw(double d , int f, int n) 
{  
  std::string s = realnum(d) ? ( f < 0 ? dbl2str(d,-f) : dbl2str_fixed(d,f) ) : "NA";
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

inline std::string sw(int i , int n) 
{
  std::string s = realnum(i) ? int2str(i) : "NA";
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

inline void NoMem()
{
  std::cerr << "*****************************************************\n"
            << "* FATAL ERROR    Exhausted system memory            *\n"
            << "*                                                   *\n"
            << "* You need a smaller dataset or a bigger computer...*\n"
            << "*                                                   *\n"
            << "* Forced exit now...                                *\n"
            << "*****************************************************\n\n";
  exit(1);
}


inline std::string itoa(int value, int base) {
  enum { kMaxDigits = 35 };
  std::string buf;
  buf.reserve( kMaxDigits ); // Pre-allocate enough space.
	
  // check that the base if valid
  if (base < 2 || base > 16) return buf;
	
  int quotient = value;
	
  // Translating number to string with base:
  do {
    buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
    quotient /= base;
  } while ( quotient );
	
  // Append the negative sign for base 10
  if ( value < 0 && base == 10) buf += '-';
	
  std::reverse( buf.begin(), buf.end() );
  return buf;	
}

inline bool realnum(double d)
{
  double zero = 0;
  if (d != d || d == 1/zero || d == -1/zero) 
    return false;
  else
    return true;
}


inline std::string int2str(int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

inline std::string longint2str(long int n)
{ 
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

inline std::string dbl2str(double n, int prc)
{
  std::ostringstream s2;
  if ( prc > 0 ) 
    s2.precision(prc);
  s2 << n;
  return s2.str();
}


inline std::string dbl2str_fixed(double n, int prc)
{
  std::ostringstream s2;
  s2 << std::setiosflags( std::ios::fixed );
  if ( prc > 0 ) 
    s2.precision(prc);
  s2 << n;
  return s2.str();
}


}
#endif
