/****************************************************************************************
 * File: clustering_quality.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 26/01/2015

 ***************************************************************************************/
#ifndef SAMOGWAS_CLUSTERING_QUALITY_HPP
#define SAMOGWAS_CLUSTERING_QUALITY_HPP

#include <set>
#include <vector>
#include <algorithm>  // std::sort

#include "partition.hpp"
namespace samogwas
{

struct ClusteringQualityCriteria {
  virtual double compute( const Clustering& clustering ) const {
    Partition p(clustering);
    return compute(p); 
  }

  virtual double compute(const Partition& part) const = 0;

  virtual double operator()(const Partition& partition) const {
    return compute(partition);
  }
  
};


struct SpatialCoherenceCriteria: public ClusteringQualityCriteria {
  virtual double compute(const Partition& partition) const {
    auto clustering = partition.to_clustering();
    return compute(partition, clustering);
  }

  virtual double compute(const Clustering& clustering) const {
    Partition partition(clustering);
    return compute(partition, clustering);
  }

  virtual double compute(const Partition& partition, const Clustering& clustering ) const {    
    double nbrChanges = nbrLabelChanges(partition);
    double expected = expectedNbrLabelChanges(clustering);
    double denominator = partition.nbrItems();

    printf("nbrItems: %d - expected: %f, nbrChanges: %f\n", partition.nbrItems(), expected, nbrChanges);
    return denominator == 0 ? 0 : 1.0 - (nbrChanges-expected)/denominator;
  }

  double nbrLabelChanges(const Partition& partition) const {
    auto nbrItems = partition.nbrItems();
    double result = 0.0;
    for (size_t i = 0; i < nbrItems-1 ;++i) {
      auto curr = partition.getLabel(i), nxt = partition.getLabel(i+1);
      if (curr != nxt) {
        result += 1;
      }
    }

    return result;
  }

  double expectedNbrLabelChanges(const Clustering& clustering) const {
    double sos = 0.0;
    double nbrItems = 0;
    for ( auto clt: clustering ) {
      sos += clt.size()*clt.size();
      nbrItems += clt.size();
      // printf("clt-size: %d\n", clt.size());
    }
    printf("expected computing: nbrItems: %f, sos: %f\n", nbrItems, sos);
    return (sos == 0.0 || nbrItems == 0) ? 0 : nbrItems - sos/nbrItems;
  }
};



} // namespace samogwas ends here. 

/****************************************************************************************/
#endif // SAMOGWAS_CLUSTERING_QUALITY_HPP
