/****************************************************************************************
 * File: rand.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 09/01/2015

 ***************************************************************************************/
#ifndef SAMOGWAS_RAND_HPP
#define SAMOGWAS_RAND_HPP

#include "compare_measure.hpp"

namespace samogwas
{

struct RandIndex: public ComparisonMeasure {
  virtual double compute( const Clustering& clusteringA, const Clustering& clusteringB ) const {
    Partition pA (clusteringA), pB(clusteringB);
    return compute( pA, pB );
  }

  //  
  virtual double compute( const Partition& pA, const Partition& pB ) const {    
    // assert( pA.nbr_objects() == pB.nbr_objects() );
    assert( pA.nbrItems() == pB.nbrItems() );

    unsigned n_11 = 0, n_00 = 0; // n_01 = 0, n_10 = 0,
    size_t N = pA.nbrItems();
    
    for ( size_t i = 0; i < N; ++i )
      for ( size_t j = i+1; j < N; ++j ) {
        int c_iA = pA.getLabel(i), c_iB = pB.getLabel(i), c_jA = pA.getLabel(j), c_jB = pB.getLabel(j);
        if ( c_iA == c_jA && c_iB == c_jB ) ++n_11;
        else if ( c_iA != c_jA && c_iB != c_jB) ++n_00;     
      }      

    return (n_11+n_00)*2 / double(N*(N-1));
  }

  virtual std::string name() const { return "RandIndex"; }

};


///////////////////////////////////////////////////////////////////////////////////////////////
struct AdjustedRandIndex: public RandIndex {
  
  virtual double compute( const Partition& pA, const Partition& pB ) const {
    assert( pA.nbrItems() == pB.nbrItems() );
    const size_t R = pA.nbrClusters(), S = pB.nbrClusters(), N = pA.nbrItems();
    // printf("%d vs %d - R: %d, S: %d, N: %d\n", pA.nbrItems(), pB.nbrItems(), R, S, N);

    // printf("1\n");     
    std::vector< std::vector<int> > conTab(R, std::vector<int>(S, 0));   
    std::vector<int> sumA(R, 0), sumB(S,0);    
    for ( size_t i = 0; i < N; ++i ) {
      int cA = pA.getLabel(i), cB = pB.getLabel(i);
      ++conTab[cA][cB];
      ++sumA[cA];
      ++sumB[cB];
    }
    // printf("2\n");     

    double expected_idx = 0, idx = 0, max_idx = 0, rs = 0;
    for ( size_t cA = 0; cA < R; ++cA ) {
      for ( size_t cB = 0; cB < S; ++cB ) {
        int val = conTab[cA][cB];
        idx += double(val*(val-1)) / 2;        
      }
    }
    // printf("3\n");     

    double expected_idxA = 0, expected_idxB;
    for ( auto& a: sumA ) {
      expected_idxA += double(a*(a-1))/2;
    }

    for ( auto& b: sumB ) {
      expected_idxB += double(b*(b-1))/2;
    }
    // printf("4\n");     

    expected_idx = 2 * expected_idxA * expected_idxB / (N*(N-1));
    max_idx = 0.5*( expected_idxA + expected_idxB );
    rs = ( max_idx == expected_idx ) ? 0.0 : ( idx - expected_idx ) / ( max_idx - expected_idx );
    return rs;
  }

  virtual std::string name() const { return "AdjustedRandIndex"; }

};

} // namespace samogwas ends here. 

/****************************************************************************************/
#endif // SAMOGWAS_RAND_HPP
