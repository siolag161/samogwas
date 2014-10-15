/****************************************************************************************
 * File: logs_utils.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 15/10/2014

 ***************************************************************************************/
#ifndef UTILITY_LOGS_UTILS_HPP
#define UTILITY_LOGS_UTILS_HPP

#include <chrono>

#include <string>
#include <typeinfo>

#ifdef __GNUG__
#include <cstdlib>
#include <memory>
#include <cxxabi.h>
#endif

namespace utility
{

#ifdef __GNUG__
static std::string demangle(const char* name) {

    int status = -4; // some arbitrary value to eliminate the compiler warning

    // enable c++11 by passing the flag -std=c++11 to g++
    std::unique_ptr<char, void(*)(void*)> res {
        abi::__cxa_demangle(name, NULL, NULL, &status),
        std::free
    };
    return (status==0) ? res.get() : name ;
}

#else
// does nothing if not g++
static std::string demangle(const char* name) {
  return name;
}
#endif

template <class T>
std::string type_name(const T& t) {
  return demangle(typeid(t).name());
}

////////////////////////////////////////////////

template<typename TimeT = std::chrono::milliseconds>
struct measure
{     
  template<typename F, typename ...Args>
  static typename TimeT::rep execution(F& func, Args&&... args)
  {
    auto start = std::chrono::system_clock::now();

    // Now call the function with all the parameters you need.
    func(std::forward<Args>(args)...);
    auto duration = std::chrono::duration_cast<TimeT>(std::chrono::system_clock::now() - start);

    return duration.count();
  }

  template<typename F, typename ...Args>
  static void log_execution( std::string name, std::string unit, F& func, Args&&... args )
  {
    std::cout << "starting performing: " << name
              << std::endl;
    auto time_taken = measure<TimeT>::execution(func, args...);
    std::cout << "end performing. took: " << time_taken << " " << unit << "(s)." << std::endl;
  }
};


} // namespace utilityends here. utility

/****************************************************************************************/
#endif // UTILITY_LOGS_UTILS_HPP
