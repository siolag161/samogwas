/****************************************************************************************
 * File: permutation_test.hpp
 * Description: This module provides a generic procedure dedicated to the correction for multiple tests.
 * ************ The distribution of the statistic under the null hypothesis is obtained by permutations.
 * 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 25/06/2014

 ***************************************************************************************/
#ifndef STATS_PERMUTATION_TEST_HPP
#define STATS_PERMUTATION_TEST_HPP

#include <boost/random.hpp> // boost::mt19937, boost::uniform_int
#include <random> // std::random_shuffle

#include <cmath>
#include <chrono> // std::chrono::system_clock
#include <algorithm> // std::min, std::max
#include <omp.h> // openMP pragmas

// #include "VecTypeype.hpp"
namespace stats
{

struct CollectionPermute {

    CollectionPermute( unsigned long seed = 1 ) {
    rng.seed(seed);
  }

  // returns the current time in nanoseconds.
  static unsigned long currentTime() {
    unsigned long time =
        std::chrono::system_clock::now().time_since_epoch() / 
        std::chrono::nanoseconds(1);
    return time;
  }
  
  /** Internal functor which returns a random unsigned integer in a given interval [0,upperLim[. 
   *  It utilizes the uniform distribution for generating the random number.
   */
  struct Rand: std::unary_function<unsigned, unsigned> {
    
    Rand(boost::mt19937 &s) : state(s) {} 

    /////////////////////////////////////////
    // returns a random unsigned integer in a given interval [0,upperLim[.
    unsigned operator()(unsigned upperLim) {
      boost::uniform_int<> rng(0, upperLim - 1);
      return rng(state);
    }
    
    boost::mt19937 &state;

  };

  ///////////////////////////////////////////////////////////////
  // Permutation with the default random number generator
  template<typename p_valueype>
  void operator()(VecType& vec) {    
    std::random_shuffle(vec.begin(), vec.end()));
  }
  
  // Permutation with a specific state 
  template<typename VecType>
  void operator()(VecType& vec, boost::mt19937 &state) {    
    Rand rand(state);
    std::random_shuffle(vec.begin(), vec.end()), rand);
  }

  // Permutation with the current state
  template<typename VecType>
  void operator()(VecType& vec) {
    Rand rand(CollectionPermute::rng); // current state
    std::random_shuffle(vec.begin(), vec.end(), rand );
  }
 private:
  boost::mt19937 rng; // random number generator of type Mersenne Twister 19937
};

//////////////////////////////////////////////////////////////////////////////////
// Returns the p-value of observing a value less than or equal to v. 
// This value is sampled from the true theoretical distribution D and dist is an empirical sample of D.
template<typename T>
double p_value( const T v, const std::vector<T>& dist ) {
  double count = 0.0;
  for (auto& val: dist ) {
    if ( val < v) ++count;
  }
  return count / dist.size();
}

//////////////////////////////////////////////////////////////////////////////////
/**  *todo: remove graph from the function */
template< class Matrix, class vector >
void performTesting( std::vector< std::vector< std::vector<double> > >& dists,
                     std::vector<std::vector<double> >& result,
                     std::vector<StatTest*> statTests,
                     const Matrix& genoMat,
                     const vector& phenotype,
                     const fltm::Graph& graph,
                     const int permu )
{
  int NP = 2;
  size_t nvars = genoMat.size(), ntests = statTests.size();
  int levels = fltm::graph_height( graph );
  
  result.resize( 2*ntests,  std::vector<double>(nvars, 0.0) );
  vector pheno = phenotype;

  
  for ( size_t test = 0; test < ntests; ++ test ) {
    for ( size_t snp = 0; snp < nvars; ++snp ) {
  
      result[2*test][snp] = statTests[test]->execute( genoMat.at(snp), pheno,
                                                      graph[snp].variable.cardinality(),
                                                      NP );
      // printf("snp: %d: %f\n", snp,t][snp]  );

    }    
  }

  if ( permu > 0 ) {
    dists.resize( ntests,
                  std::vector< std::vector<double> >( (levels+1),
                                                      std::vector<double>( permu, 2.0 )));
    

    CollectionPermutate permute;
    #pragma omp parallel for
    for ( int p = 0; p < permu; ++p ) {
       #pragma omp critical
       permute(pheno);
       for ( size_t test = 0; test < ntests; ++ test ) {
         if (statTests[test]->name == "Fisher") omp_set_num_threads(1);
         for ( size_t snp = 0; snp < nvars; ++snp ) {
           double pval = statTests[test]->execute( genoMat.at(snp), pheno, graph[snp].variable.cardinality(), NP);
           int level = graph[snp].level;          
           #pragma omp critical
           dists[test][level][p]= std::min( dists[test][level][p], pval );
         }    
      }       
    }

    for ( size_t snp = 0; snp < nvars; ++snp ) {
      int level = graph[snp].level;
      for ( size_t test = 0; test < ntests; ++ test ) {
         result[2*test+1][snp] = p_value( result[2*test][snp], dists[test][level] );
       }   
     }
  }
}


/////////////////////////////////
/** todo: change to list (with name and such)
 *
 //  */
template< class Matrix, class VecType >
void performTesting( std::vector<double>& dist,
                     std::vector<double>& result,
                     StatTest* statTest,
                     const std::vector<int> candidates,
                     const Matrix& genoMat,
                     const VecType& phenotype,
                     const fltm::Graph& graph,
                     const int permu )
{
  int NP = 2;
  size_t nvars = candidates.size();
  result.resize( 2*nvars, 0.0 );
  vector pheno = phenotype;
  for ( size_t snp = 0; snp < nvars; ++snp ) {
    int candidate = candidates.at(snp);
    result[2*snp] = statTest->execute( genoMat.at(candidate), pheno,
                                       graph[candidate].variable.cardinality(),
                                       NP );
  }    
  

  if ( permu > 0 ) {
    dist.resize(permu, 2.0);
    CollectionPermutate permute;
    #pragma omp parallel for  
    for ( int p = 0; p < permu; ++p ) {
       #pragma omp critical
       permute(pheno);
       if (statTest->name == "Fisher") omp_set_num_threads(1);
       for ( size_t snp = 0; snp < nvars; ++snp ) {
         int candidate = candidates.at(snp);
         double pval = statTest->execute( genoMat.at(candidate), pheno, graph[candidate].variable.cardinality(), NP);
         #pragma omp critical
         dist[p]= std::min( dist[p], pval );
       }               
    }

    for ( int snp = 0; snp < nvars; ++snp ) {
      int candidate = candidates.at(snp);
      result[2*snp+1] = p_value( result[2*snp], dist );       
    }
  }
}


} // namespace stats ends here. 


/****************************************************************************************/
#endif //STATS_PERMUTATION_TEST_HPP
