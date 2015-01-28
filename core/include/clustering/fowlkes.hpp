/****************************************************************************************
 * File: folkwes.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 09/01/2015

 ***************************************************************************************/
#ifndef SAMOGWAS_FOLKWES_HPP
#define SAMOGWAS_FOLKWES_HPP

#include "compare_measure.hpp"
namespace samogwas
{

struct FowlkesIndex: public ComparisonMeasure {
  virtual double compute( const Clustering& clusteringA, const Clustering& clusteringB ) const {
    Partition pA (clusteringA), pB(clusteringB);
    return compute( pA, pB );
  }

  //  
  virtual double compute( const Partition& pA, const Partition& pB ) const {
    // printf("%d vs %d\n", pA.nbrItems(), pB.nbrItems());
    assert( pA.nbrItems() == pB.nbrItems() );
    int n_11 = 0, n_00 = 0, n_01 = 0, n_10 = 0;
    size_t N = pA.nbrItems();
    
    for ( size_t i = 0; i < N; ++i )
      for ( size_t j = i+1; j < N; ++j ) {
        int c_iA = pA.getLabel(i), c_iB = pB.getLabel(i), c_jA = pA.getLabel(j), c_jB = pB.getLabel(j);
        if ( c_iA == c_jA && c_iB == c_jB ) ++n_11;
        else if ( c_iA != c_jA && c_iB != c_jB) ++n_00;
        else if ( c_iA == c_jA && c_iB != c_jB) ++n_01;
        else if ( c_iA != c_jA && c_iB == c_jB) ++n_10;
      }      

    return std::sqrt( double(n_11*n_11)/ ((n_11+n_01)*(n_11+n_10)) );
  }

  virtual std::string name() const { return "Fowlkes-Mallows-Index"; }
};


} // namespace samogwas ends here. 

/****************************************************************************************/
#endif // SAMOGWAS_FOLKWES_HPP
