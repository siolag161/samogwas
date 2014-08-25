/****************************************************************************************
 * File: em_helper.hpp
 * Description: Some helpers for Customized EM Algorihtms
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 05/08/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_EM_HELPER_HPP
#define SAMOGWAS_EM_HELPER_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

namespace samogwas
{

typedef std::vector<double> Distribution;
typedef std::vector< std::vector<double> > ProbTable;  
typedef std::vector< std::vector< std::vector<double> > > CondProbTable;

typedef std::vector< std::vector<int> > Matrix;  

struct NV_EM_log_likelihood {
  
  double likelihood( const Matrix& data,
                     const unsigned i,
                     const ProbTable& theta,
                     const Distribution& pY,
                     const CondProbTable& pXY ) const;
  
  double operator()( const Matrix& data,
                     const ProbTable& theta,
                     const Distribution& pX,
                     const CondProbTable& pXY ) const;  
};

///////////////////////////////////////////////////////////

/** Takes a collection of values and normalize them to get a distribution ( scaling to [0-1] and sum to 1.0 )
 *
 */
struct Density_Normalize {
  template< typename RealCollection >
  void operator()( RealCollection& rc ) {
    double sum = 0.0; for ( auto& item: rc ) sum += item;
    for ( int i = 0; i < rc.size(); ++i ) rc[i] /= sum;
  }
};


/** This provides a convenient way to computes log(X) by taking care of the case where X = 0
 *
 */
template<typename T>
double m_log( const T t) {
  return t <= 0 ? 0.0 : std::log(t);
}

/** This provides a convenient way to generate a random integer
 *
 */
struct Int_Rand {
  unsigned operator()(unsigned i) {
    std::uniform_int_distribution<> dist(0, i-1);
    return dist(gen);
  }

  Int_Rand(): gen(std::random_device()()) {}  
  std::mt19937 gen;
};

struct Rand_Density {
  template<typename T>
  void operator()( std::vector<T>& prob, int sz, int max_val ) {
    Int_Rand int_rand;
    prob.resize(sz, 1.0);
    for ( int i = 0; i < sz; ++i ) {
      prob[i] = int_rand(max_val);
    }

    Density_Normalize()(prob);
  }
};

//////////////////////////////////////////

/** This provides a simple implementation for Kullback–Leibler divergence for comparing 
 *  a distribution to its estimtation
 */
template< typename PC, typename QC >
double KL( const PC& P, const QC& Q ) {
  assert( P.size() == Q.size() );
  double rs = 0.0;
  for ( int i = 0; i < P.size(); ++i ) {
    rs += Q[i] == 0 ? 0.0 : m_log(P[i]/Q[i])*P[i];
  }
  return rs;
}


} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_EM_HELPER_HPP
