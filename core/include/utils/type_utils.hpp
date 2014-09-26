/*********************************************************************
 * File: type_utils.hpp
 * Description:  This module contains several helper methods related to type manipulation
 * ------------  and casting between different types (hence its name).
 * 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet 
 * @date: 18/10/2013
 *********************************************************************/

#ifndef UTILS_TYPE_UTILS_HPP
#define UTILS_TYPE_UTILS_HPP

#include <sstream> // stringstream
#include <string> // basic_string

namespace utility {    

/******************************************* DECLARATION *****************************************/
/**
 * Constructs enum from index.
 */
// template<typename EnumType>
// EnumType EnumOfIndex(const int& i); 

/**
 * Integer to Type trick (from Modern C++ Design)
 */
template<int I>
struct Int2Type
{
  enum {value = I};
};

      
///////////////////////// generic from string and to string ///////////////////////////////
/**
 * Generic functor to convert the string (normal, wide string)... to ConvertType
 */
template<typename CharType, typename ConvertType> 
struct fromString: public std::unary_function<CharType, ConvertType>
{
  ConvertType operator()(const std::basic_string<CharType>& str);
 private:
  std::stringstream ss;
};


/**
 * Handles the conversion from string to string. 
 */
template<typename CharType> 
struct fromString<CharType, std::string >
    : public std::unary_function<CharType, std::string >
{
  std::string operator()(const std::string& str);
};


/**
 * Generic functor to convert ConvertType to string (normal, wide string).
 */
template<typename CharType, typename ConvertType> 
struct toString: public std::unary_function<ConvertType, std::basic_string<CharType> >
{
  std::basic_string<CharType> operator()(const ConvertType&);
  
 private:
  std::stringstream ss;
};


/**
 * Handles the conversion from string to string. 
 */
template<typename CharType> 
struct toString<CharType, std::string >
    : std::unary_function<CharType, std::string >
{
  std::string operator()(const std::string&);

};
}

/************************************************** IMPLEMENTATION **************************************************/

namespace utility {
   
template<typename EnumType>
EnumType EnumOfIndex(const int& i)
{
  return static_cast<EnumType>(i);
}

////////////////////// FROM STRING AND TO STRING //////////////////////

///////////////////////////////////////////////////////////////////////////
/**
 */
template<typename CharT, typename ConvertT> 
ConvertT fromString<CharT, ConvertT>::operator()(const std::basic_string<CharT>& str)
{
  ConvertT result(0);

  if (!(ss<<str))
  {
    //throw std::runtime_error("cannot convert")
    ss.clear();
  }
  if (!(ss>>result)) 
  {
    //throw std::runtime_error("cannot convert")
    ss.clear();
  }
  ss.clear();
  return result;
}

/** Specific case for fromString (string to string).
 */
template<typename CharT> 
std::string fromString<CharT, std::string >::operator()
    (const std::string& str)
{
  return str;
} 

//////////////////////////////////////////////////////////////////
template<typename CharT, typename ConvertT> 
std::basic_string<CharT> toString<CharT, ConvertT>::operator()(const ConvertT& t) {
  std::basic_string<CharT> result;

  if (!(ss<<t)) {
    //throw std::runtime_error("cannot convert")
    ss.clear();
  }
  if (!(ss>>result)) 
  {
    //throw std::runtime_error("cannot convert")
    ss.clear();
  }
  ss.clear();
  return result;
}


/** Specific case for toString (string to string).
 */
template<typename CharT> 
std::string toString<CharT, std::string >::operator()
    (const std::string& str)
{
  return str;
}

}

#endif // UTILS_TYPE_UTILS_HPP
