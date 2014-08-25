
#include "em/em_helper.hpp"
#include "utils/matrix_utils.hpp"
#include <cmath>

#include <cstdio>
namespace samogwas
{

double NV_EM_log_likelihood::likelihood( const Matrix& data,
                                         const unsigned i,
                                         const ProbTable& theta,
                                         const Distribution& pY,
                                         const CondProbTable& pXY) const {
  double llh = 0.0;
  const std::vector<int>& X = data[i];
  unsigned K = pY.size(), P = X.size();

  for ( unsigned y = 0; y < K; ++y ) {
    double llh_y = m_log(pY[y]);
    for ( int p = 0; p < P; ++p ) {
      int x = X[p];
      llh_y += m_log( pXY[y][p][x] );
    }
    llh_y *= theta[i][y];
  }
  
  return llh;
}

///////////////////////////////////////////////////////////////////

double NV_EM_log_likelihood::operator()( const Matrix& data,
                                         const ProbTable& theta,
                                         const Distribution& pY,
                                         const CondProbTable& pXY ) const {
  double result = 0.0;
  unsigned N = utility::nrows(data);

  for ( unsigned i = 0; i < N; ++i ) {
    const std::vector<int>& X = data[i];
    unsigned K = pY.size(), P = X.size();

    for ( unsigned y = 0; y < K; ++y ) {
      double llh_y = m_log(pY[y]);
      for ( int p = 0; p < P; ++p ) {
        int x = X[p];
        if (pXY[y][p][x]) {
          llh_y += m_log( pXY[y][p][x] );
        }
      }
      llh_y *= theta[i][y];
      result += llh_y;      
    }    
          
  }
  // printf("result: %f\n", result );
  return result;
}






} // namespace samogwasends here. samogwas
