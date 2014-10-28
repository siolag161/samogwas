/**                     NOT REVIEWED
 *
 * ****************************************************************************
 * File: multi_em.cpp
 * Description:
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 03/08/2014

 ***************************************************************************************/

#define NOT_MISSING 1
#define MISSING 0
#define DEFAULT_VALUE = -1

#include "utils/matrix_utils.hpp"
#include "em/multi_em.hpp"
#include "em/em_helper.hpp"

namespace samogwas
{

void NV_EM::clear_values() {
  theta.clear(); theta = std::vector< std::vector<double> >();
  pY.clear(); pY = std::vector<double>();
  pYX.clear(); pYX = std::vector< std::vector< std::vector<double> > >();
}

///////////////////////////////////////////////////////////////
void NV_EM::init_values( int N, int K, const Variables& variables ) {
  Rand_Density rand_density;
  rand_density( pY, K, N);
  for ( int i = 0; i < N; ++i ) {
    std::vector<double> probTab;
    rand_density( probTab, K, N);
    theta.push_back( std::vector<double>( K, 1.0/K ) );
  }

  unsigned int P = variables.size();
  for ( int i = 0; i < K; ++i ) {
    std::vector< std::vector<double> > probTables;
    probTables.reserve(P);
    for ( int i = 0; i < P; ++i ) {
      unsigned int card = variables[i].cardinality();
      std::vector<double> probTab;
      rand_density( probTab, card, N);
      probTables.push_back(probTab);
    }
    pYX.push_back(probTables);
  }
}


void NV_EM::run( ResultEM& result,
                 const Variable& latentVar,
                 const Variables& variables,
                 const NV_EM::Matrix& data,
                 const double threshold,
                 const std::vector< std::vector<bool> > & defTable ) {
  
  int N = utility::nrows(data), K = latentVar.cardinality();
  init_values( N, K, variables );
  NV_EM_log_likelihood log_likelihood;
  
  double prev_llh = 0.0;
  double curr_llh = log_likelihood(data, theta, pY, pYX);

  int step = 0;
  // print_tabs( theta, pY, pYX );
  while ( ++step < 200  ) {
    step_E( theta, pY, pYX, data );
    step_M( pY, pYX, theta, data );

    prev_llh = curr_llh;
    curr_llh = log_likelihood(data, theta, pY, pYX);

    if ( std::abs( prev_llh - curr_llh) < threshold ) break;
  }
  std::cout << "takes: " << step << " steps\n";
  // print_tabs( theta, pY, pYX );
  
}

/**
 *
 */
void NV_EM::imputeLatent( ResultEM& result,                 
                          const plSymbol& latentVar,
                          const NV_EM::Matrix& dataTable,
                          EMLearner& bestModel,
                          plMatrixDataDescriptor<int>& dataDesc ) {
  
}

/**
 *
 */
void NV_EM::step_E( std::vector< std::vector<double> >& theta,
                    const std::vector<double>& pY,
                    const std::vector< std::vector< std::vector<double> > >& pYX,
                    const NV_EM::Matrix& data ) {
  int N = utility::nrows(data);
  for ( int i = 0; i < N; ++i ) {
    double sum = 0.0;
    for ( int y = 0; y < pY.size(); ++y ) {
      theta[i][y] = m_log(pY[y]);
      for ( int p = 0; p < pYX[y].size(); ++p ) {
        int x = data[i][p];
        theta[i][y] += m_log(pYX[y][p][x]);
      }
      theta[i][y] = std::exp(theta[i][y]);
      sum += theta[i][y];
    }
    for ( int y = 0; y < pY.size(); ++y ) {
      theta[i][y] /= sum;
    }
  }
}

/**
 *
 */
void NV_EM::step_M( std::vector<double>& pY,
                    std::vector< std::vector< std::vector<double> > >& pYX,
                    const std::vector< std::vector<double> >& theta,
                    const NV_EM::Matrix& data ) {
  int N = theta.size();
  for ( int y = 0; y < pY.size(); ++y ) { // priors
    pY[y] = 0.0;
    for ( int i = 0; i < N; ++i) {
      pY[y] += theta[i][y];
    }
    pY[y] /= N;
  }

  for ( int y = 0; y < pY.size(); ++y ) {
    for ( int p = 0; p < pYX[y].size(); ++p ) {
      double sum = 0.0;
      std::vector<double> xvals( pYX[y][p].size(), 0.0 );
      std::vector<int> counts( pYX[y][p].size(), 0 );

      for ( int i = 0; i < N; ++i ) {
        int x = data[i][p];
        xvals[x] += theta[i][y];
        sum += theta[i][y];
        counts[x]++;
      }      
    } // for p 
  } // for y

}

} // namespace samogwasends here. samogwas
