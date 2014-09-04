/****************************************************************************************
 * File: matrix_utils.hpp
 * Description: This modules regroups several commonly used mathods for working with
 * ------------ matrix structure ( a vector-of-vector like objects ).
 * ------------ We have methods dealing with cardinality and computing like Tranpose
 * ------------ a matrix are included.
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 29/04/2014

 ***************************************************************************************/
#ifndef UTILS_MATRIX_HPP
#define UTILS_MATRIX_HPP

#include <vector>
namespace utility // common namespace in this sub-repo
{

/** Return number of rows of a matrix-like data structure. We suppose the input matrix 
 *  has a row-major skeleton ( hence the need for the size() method).
 */
template<typename MatrixT>
unsigned nrows(const MatrixT& mat) {
  return mat.size();
}

/** Return number of column of a matrix-like data structure. Same hypothesis as above.
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

} // namespace Utilityends here. Utility



/****************************************************************************************/
#endif // UTILS_MATRIX_HPP
