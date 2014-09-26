/****************************************************************************************
 * File: csv_parser.hpp
 * Description: This module is responsible for reading a CSV-like file and outputing the result into 
 * ------------ a Matrix-type object. Since the method is templated, a row in the input file can contain mixed types:
 * ------------ for example: Tom; 30.32; Male. All the rows have the same structure.
 * 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet @TODO: modify this part
 * @date: 25/03/2014

 ***************************************************************************************/
#ifndef UTILS_CSV_PARSER_HPP
#define UTILS_CSV_PARSER_HPP

#include <iterator> // std::input_iterator_tag 
#include <fstream> // std::getline
#include <sstream> // std::lineStream
#include <iostream> // std::istream
#include <vector>
#include <string>

namespace utility // the common namesapce for utils files @todo: utility -> utils
{

/** If in the input file some fields have different types, T = std::string must be used.
 *  Otherwise, T is specified as the common type of all the fields. 
 */
template<class T>
class CSVRow
{
 public:
  // retrieves the stored data at the index position of the row.
  T const& operator[](std::size_t index) const {
    return m_data[index];
  }

  T const at(std::size_t index) const
  {
    return m_data[index];
  }
  
  std::size_t size() const
  {
    return m_data.size();
  }

  /** calls this function to read the next row (if any) and update the current stored data with this row. 
   */
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
      if ( lineStream.peek() == ',') // skips the separator character // @todo: ',' -> parameter
        lineStream.ignore();
    }
  }
  
 private:
  std::vector<T> m_data;
};

template< typename T>
std::istream& operator>>(std::istream& str, CSVRow<T>& data)
{
  data.readNextRow(str);
  return str;
}


/** The iterator itself. It implements all the common iterator interface methods
 *  std::ifstream infile("in.txt");"// it.txt contains 2 lines: 1,2 and 3,4
 *  CSVIterator<std::string> it(infile);
 *  *(++it) returns 3,4 and now points to 3,4 while *(it++) returns 1,2 and also points to 3,4<
 */
template<typename T>
class CSVIterator
{    
 public:
  /** Structures used by iterator interface
   *
   */
  typedef std::input_iterator_tag iterator_category;
  typedef CSVRow<T> value_type;
  typedef std::size_t difference_type;
  typedef CSVRow<T>* pointer;
  typedef CSVRow<T>& reference;

  utility::CSVIterator<T> matrixLine(matrixFile);
  
  for( ; matrixLine != utility::CSVIterator<T>(); ++matrixLine ) {
  utility::CSVIterator<T> matrixLine(matrixFile);
  
  for( ; matrixLine != utility::CSVIterator<T>(); ++matrixLine ) {
  CSVIterator(std::istream& str): m_str(str.good() ? &str : NULL) { ++(*this); }
  CSVIterator(): m_str(NULL) {}

  // Pre Increment 
  CSVIterator& operator++() { if (m_str) { (*m_str) >> m_row;m_str = m_str->good()?m_str:NULL;}return *this;}

  // Post increment (copy)
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

} // namespace utility ends here.

/****************************** IMLEMENTATION BELOW THIS POINT **************************/
namespace utility
{

/** This template specialization is meant to override the default template implementation above (for readNextRow)
 *  when the template parameter T is std::string.
 */
template<>
class CSVRow<std::string>
{
 public:
  std::string const& operator[](std::size_t index) const
  {
    return m_data[index];
  }

  std::string const at(std::size_t index) const
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
