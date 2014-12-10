/****************************************************************************************
 * file: permutation_test.hpp
 * description: this module provides a generic procedure dedicated to the correction for multiple tests.
 * ************ the distribution of the statistic under the null hypothesis is obtained by permutations.
 *
 * @author: duc-thanh phan siolag161 (thanh.phan@outlook.com), under the supervision of christine sinoquet
 * @date: 25/06/2014

 ***************************************************************************************/
#ifndef stats_permutation_test_hpp
#define stats_permutation_test_hpp

#include <boost/random.hpp> // boost::mt19937, boost::uniform_int
#include <random> // std::random_shuffle

#include <cmath>
#include <chrono> // std::chrono::system_clock
#include <algorithm> // std::min, std::max
#include <omp.h> // openmp pragmas
#include "statistics/association_test.hpp"
// #include "vectypeype.hpp"
namespace stats
{

struct CollectionPermute {

  CollectionPermute( unsigned long seed = 1 ) {
    rng.seed(seed);
  }

  // returns the current time in nanoseconds.
  static unsigned long currenttime() {
    unsigned long time =
        std::chrono::system_clock::now().time_since_epoch() /
        std::chrono::nanoseconds(1);
    return time;
  }

  /** internal functor which returns a random unsigned integer in a given interval [0,upperlim[.
   *  it utilizes the uniform distribution for generating the random number.
   */
  struct rand: std::unary_function<unsigned, unsigned> {

    rand(boost::mt19937 &s) : state(s) {}

    /////////////////////////////////////////
    // returns a random unsigned integer in a given interval [0,upperlim[.
    unsigned operator()(unsigned upperLim) {
      boost::uniform_int<> rng(0, upperLim - 1);
      return rng(state);
    }

    boost::mt19937 &state;

  };

  ///////////////////////////////////////////////////////////////
  // permutation with the default random number generator
  template<typename VecType>
  void operator()(VecType& vec) {
    std::random_shuffle(vec.begin(), vec.end());
  }

  // permutation with a specific state
  template<typename VecType>
  void operator()(VecType& vec, boost::mt19937 &state) {
    rand rand(state);
    std::random_shuffle(vec.begin(), vec.end(), rand);
  }

  // permutation with the current state
  template<typename VecType>
  void operator()(VecType vec) {
    rand rand(rng); // current state
    std::random_shuffle(vec.begin(), vec.end(), rand );
  }
 private:
  boost::mt19937 rng; // random number generator of type mersenne twister 19937
};

//////////////////////////////////////////////////////////////////////////////////
// returns the p-value of observing a value less than or equal to v. 
// this value is sampled from the true theoretical distribution d and dist is an empirical sample of d.
template<typename T>
double p_value( const T v, const std::vector<T>& dist ) {
  double count = 0.0;
  for (auto& val: dist ) {
    if ( val < v) ++count;
  }
  return count / dist.size();
}
//////////////////////////////////////////////////////////////////////////////////
template< class Matrix, class Vector >
void permutationProcedure( std::vector<double> &distri,
                           std::vector<double> &pvals,
                           StatTest* statTest,
                           const Matrix &xData,
                           const Vector &yData,
                           const std::vector<size_t> &xCandidates,
                           const std::vector<int> &xCardinalities,
                           const int yCardinality,
                           const int permu = 2) {
  assert(xCardinalities.size() == xCandidates.size());

  size_t nvars = xCandidates.size();
  pvals.resize(nvars, 0.0);

  for ( size_t c = 0; c < nvars; ++c) {
    size_t var = xCandidates[c];
    pvals[c] = statTest->execute( xData.at(var), yData,
                                  xCardinalities[c], permu );
  }

  if ( permu > 0 ) {
    Vector localYData = yData;
    distri.resize(permu, 2.0);
    for ( int p = 0; p < permu; ++p ) {
      permute(localYData);
      for ( size_t c = 0; c < nvars; ++c) {
        size_t var = xCandidates[c];
        double pval = statTest->execute( xData.at(var), localYData,
                                         xCardinalities[c], permu );
        distri[p] = std::min( distri[p], pval );
      }
                
    }

    for ( size_t c = 0; c < nvars; ++c) {
      size_t var = xCandidates[c];
      pvals[c]= p_value( pvals[c], distri );

    }
  }
}
//////////////////////////////////////////////////////////////////////////////////
/**  *todo: remove graph from the function *//*

     template< class Matrix, class vector >
     void performtesting( std::vector< std::vector< std::vector<double> > >& dists,
     std::vector<std::vector<double> >& result,
     std::vector<StatTest*> StatTests,
     const Matrix& genomat,
     const vector& phenotype,
     const samogwas::graph& graph,
     const int permu )
     {
     int np = 2;
     size_t nvars = genomat.size(), ntests = StatTests.size();
     int levels = fltm::graph_height( graph );

     result.resize( 2*ntests,  std::vector<double>(nvars, 0.0) );
     vector pheno = phenotype;


     for ( size_t test = 0; test < ntests; ++ test ) {
     for ( size_t snp = 0; snp < nvars; ++snp ) {

     result[2*test][snp] = statTest->execute( genomat.at(snp), pheno,
     graph[snp].variable.cardinality(),
     np );
     }
     }

     if ( permu > 0 ) {
     dists.resize( ntests,
     std::vector< std::vector<double> >( (levels+1),
     std::vector<double>( permu, 2.0 )));


     collectionpermutate permute;
                                             */
/*#pragma omp parallel for*//*

  for ( int p = 0; p < permu; ++p ) {
                            */
/*#pragma omp critical*//*

  permute(pheno);
  for ( size_t test = 0; test < ntests; ++ test ) {
  if (statTest->name == "fisher") omp_set_num_threads(1);
  for ( size_t snp = 0; snp < nvars; ++snp ) {
  double pval = statTest->execute( genomat.at(snp), pheno, graph[snp].variable.cardinality(), np);
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
  */
/** todo: change to list (with name and such)
 *
 //  *//*

     template< class Matrix, class vectype >
     void performtesting( std::vector<double>& dist,
     std::vector<double>& result,
     StatTest* StatTest,
     const std::vector<int> candidates,
     const Matrix& genomat,
     const vectype& phenotype,
     const fltm::graph& graph,
     const int permu )
     {
     int np = 2;
     size_t nvars = candidates.size();
     result.resize( 2*nvars, 0.0 );
     vector pheno = phenotype;
     for ( size_t snp = 0; snp < nvars; ++snp ) {
     int candidate = candidates.at(snp);
     result[2*snp] = StatTest->execute( genomat.at(candidate), pheno,
     graph[candidate].variable.cardinality(),
     np );
     }


     if ( permu > 0 ) {
     dist.resize(permu, 2.0);
     collectionpermutate permute;
     #pragma omp parallel for
     for ( int p = 0; p < permu; ++p ) {
     #pragma omp critical
     permute(pheno);
     if (StatTest->name == "fisher") omp_set_num_threads(1);
     for ( size_t snp = 0; snp < nvars; ++snp ) {
     int candidate = candidates.at(snp);
     double pval = StatTest->execute( genomat.at(candidate), pheno, graph[candidate].variable.cardinality(), np);
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
       */


} // namespace stats ends here. 


/****************************************************************************************/
#endif //stats_permutation_test_hpp
