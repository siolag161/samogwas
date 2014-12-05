/****************************************************************************************
 * File: property.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 07/11/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_PROPERTY_HPP
#define SAMOGWAS_PROPERTY_HPP

namespace samogwas
{

/** @todo: description
 */
struct AbstractSettings
{  
  virtual void load(const std::string& filename) { load(filename.c_str()); }
  virtual void save(const std::string& filename) { save(filename.c_str()); }
  
  virtual void load(const char* filename) = 0;
  virtual void save(const char* filename) = 0;

  template<typename T>
  AbstractSettings& setProperty(const char* property, const T& val);

  template<typename T>
  AbstractSettings& setProperty(const std::string& property, const T& val);
};

} // namespace samogwas ends here

/****************************************************************************************/
#endif // SAMOGWAS_PROPERTY_HPP
