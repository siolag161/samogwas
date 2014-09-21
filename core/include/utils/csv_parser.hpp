/****************************************************************************************
 * File: csv_parser.hpp
 * Description: This module is responsible for reading a CSV-like files and outputs the result into 
 * ------------ a Matrix-type object. Since it's templated, the input files can have mixed types
 * ------------ for example: Tom; 30.32; Male
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet @TODO: modify this part
 * @date: 25/03/2014

 ***************************************************************************************/
#ifndef UTILS_CSV_PARSER_HPP
#define UTILS_CSV_PARSER_HPP

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace utility // the common namesapce for utils files
{

template<class T>
class CSVRow
{
 public:
  inline T const& operator[](std::size_t index) const {
    return m_data[index];
  }

  inline T const at(std::size_t index) const
  {
    return m_data[index];
  }
  
  std::size_t size() const
  {
    return m_data.size();
  }

  inline void readNextRow(std::istream& str)
  {
    std::string line;
    std::getline(str,line);

    std::stringstream lineStream(line);
    T cell;

    m_data.clear();
    while (lineStream >> cell)
    {
      m_data.push_back(cell);
      if ( lineStream.peek() == ',')
        lineStream.ignore();
    }
    
  }
 private:
  std::vector<T> m_data;
};

template< typename T>
inline std::istream& operator>>(std::istream& str, CSVRow<T>& data)
{
  data.readNextRow(str);
  return str;
}

/** The iterator itself. It implements all the common interator-interface methods
 *  in order to be used as is.
 */
template<typename T>
class CSVIterator
{    
 public:
  /** nameing convention for iterator-like common structures ( value_type, pointer, etc.)
   *
   */
  typedef std::input_iterator_tag iterator_category;
  typedef CSVRow<T> value_type;
  typedef std::size_t difference_type;
  typedef CSVRow<T>* pointer;
  typedef CSVRow<T>& reference;

  CSVIterator(std::istream& str): m_str(str.good()?&str:NULL) { ++(*this); }
  CSVIterator(): m_str(NULL) {}

  // Pre Increment
  CSVIterator& operator++() { if (m_str) { (*m_str) >> m_row;m_str = m_str->good()?m_str:NULL;}return *this;}

  // Post increment

  
  CSVIterator operator++(int) { CSVIterator tmp(*this);++(*this);return tmp;}

  // Dereferences
  CSVRow<T> const& operator*() const { return m_row;}

  // References
  CSVRow<T> const* operator->() const  { return &m_row;}

  // Equality
  bool operator==(CSVIterator const& rhs) { return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
  bool operator!=(CSVIterator const& rhs) { return !((*this) == rhs);}
  
 private:
  std::istream* m_str;
  CSVRow<T> m_row; // current row
};

} // namespace graphends here. graph

/****************************** IMLEMENTATION BELOW THIS POINT **************************/
namespace utility
{

/** CSV-Row is a interate-like ( meaning we can use it in an interation manner )
 *  which reads a row in a files and stores in a std::vector of strings
 */
template<>
class CSVRow<std::string>
{
 public:
  // retrieves the stored data at the index position of the row
  inline std::string const& operator[](std::size_t index) const
  {
    return m_data[index];
  }

  inline std::string const at(std::size_t index) const
  {
    return m_data[index];
  }
  
  std::size_t size() const
  {
    return m_data.size();
  }

  /** calls this function to clear the current row data and proceeds to read the next
   *  row. 
   */
  inline void readNextRow(std::istream& str)
  {
    std::string line;
    std::getline(str,line);

    std::stringstream lineStream(line);
    std::string cell;

    m_data.clear();
    while(std::getline(lineStream, cell, ','))
    {
      m_data.push_back(cell);
    }
    
  }
  
 private:
  // data holder
  std::vector<std::string> m_data;
};


} // namespace utils ends here.

/****************************************************************************************/
#endif // UTILS_CSV_PARSER_HPP
