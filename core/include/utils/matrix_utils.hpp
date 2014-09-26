/****************************************************************************************
 * File: matrix_utils.hpp
 * Description: This modules regroups several commonly used methods for working with a
 * ------------ matrix structure (a vector of vector-like objects).
 * 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 29/04/2014

 ***************************************************************************************/
#ifndef UTILS_MATRIX_HPP
#define UTILS_MATRIX_HPP

#include <vector>

namespace utility 
{

/** Returns the number of rows of a matrix-like data structure. We suppose the input matrix 
 *  has a row-major skeleton (hence the need for the size() method).
 */
template<typename MatrixT>
unsigned nrows(const MatrixT& mat) {
  return mat.size();
}

/** Returns the number of columns of a matrix-like data structure. Same assumption as above.
 */
template<typename MatrixT>
unsigned ncols(const MatrixT& mat) {
  return mat.empty() ? 0 : mat.at(0).size();
}

/** Takes a matrix and returns its transposed matrix.
 */
template<class T>
std::vector< std::vector<T> > Transpose( const std::vector< std::vector<T> >& mat ) {
  std::vector< std::vector<T> > result( ncols(mat), std::vector<T>( nrows(mat), 0) );

  for (unsigned row = 0; row < nrows(mat); row++) {    
    for (unsigned col = 0; col < ncols(mat); col++)
    {
      result[col][row] = mat.at(row).at(col);
    }
  }  
  return result;
}

} // namespace utility ends here.



/****************************************************************************************/
#endif // UTILS_MATRIX_HPP
