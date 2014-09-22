/****************************************************************************************
 * File: em_helper.hpp
 * Description: Some helpers // CS What is a helper? for Customized EM Algorihtms
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 05/08/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_EM_HELPER_HPP
#define SAMOGWAS_EM_HELPER_HPP

#include <vector>
#include <cmath> // CS for which function? log
#include <algorithm> // idem
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
  
  // CS no comment?
  double operator()( const Matrix& data,
                     const ProbTable& theta,
                     const Distribution& pX,
                     const CondProbTable& pXY ) const;  
};

///////////////////////////////////////////////////////////

/** Takes a collection of values and normalize<CS> them to get<obtain> a distribution ( scaling<CS> to [0-1] and sum<CS> to 1.0 )
 *
 */ // CS Is it really used?
struct Density_Normalize { // CS Why an underscore in the function name?
  template< typename RealCollection >
  void operator()( RealCollection& rc ) {
    double sum = 0.0; for ( auto& item: rc ) sum += item;
    for ( int i = 0; i < rc.size(); ++i ) rc[i] /= sum;
  }
};


/** <CS> Heterogeneous style. This provides Provides a convenient way to computes<CS> log(X) by taking care of the case where<CS when> X = 0
 *
 */
template<typename T>
double m_log( const T t) {
  return t <= 0 ? 0.0 : std::log(t);
}

/** <CS> Heterogeneous style. This provides a convenient way to generate a random integer
 *
 */
struct Int_Rand {
  unsigned operator()(unsigned i) {
    std::uniform_int_distribution<> dist(0, i-1); <CS I do not understand the uniform_int_distribution<CS> >
    return dist(gen); // CS explain
  }

  Int_Rand(): gen(std::random_device()()) {}  // CS When is it used?
  std::mt19937 gen;
};

struct Rand_Density {
  template<typename T>
  void operator()( std::vector<T>& prob, int sz, int max_val ) { // CS prob or distribution?
    Int_Rand int_rand;
    prob.resize(sz, 1.0); // CS extends the size of prob up to sz; the new elements are initialized as 1.0.
                          // CS Is it useful to care about the initialization?
    for ( int i = 0; i < sz; ++i ) {
      prob[i] = int_rand(max_val);
    }

    Density_Normalize()(prob);// CS the functiun call is not easy to understand.
  }
};

//////////////////////////////////////////

/** This provides a simple implementation for Kullback–Leibler divergence for comparing 
 *  a distribution to its <CS>estimtation
 */
 // CS formula lacking
 // D_KL(P|Q) = sum_i P(i) log (P(i)/Q(i)) 
template< typename PC, typename QC >
double KL( const PC& P, const QC& Q ) {
  assert( P.size() == Q.size() );
  double res = 0.0;
  for ( int i = 0; i < P.size(); ++i ) {
    res += Q[i] == 0 ? 0.0 : P[i] * m_log(P[i]/Q[i]);
  }
  return res;
}


} // namespace samogwasends here. samogwas 

/****************************************************************************************/
#endif // SAMOGWAS_EM_HELPER_HPP
