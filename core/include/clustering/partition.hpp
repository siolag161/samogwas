/****************************************************************************************
 * File: partition.hpp 
 * Description: A data structure (Partition) representing a partition of the data set - i.e a division of the data set into
 * -----------  non-overlapping subsets.
 * -----------  This data structure is used to store the result produced by an instanciation of 
 * -----------  AlgoClusteringInterface (algo_clustering.hpp).
 * -----------  Globally, this file provides the toolbox for algo_clustering.hpp.
 * @author: Duc-Thanh Phan Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 09/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_PARTITION_HPP 
#define SAMOGWAS_PARTITION_HPP

#include "distance/comparable.hpp"

#include <set>
#include <vector>
#include <map>

#include <iostream>

namespace samogwas
{

/** An Index depicts an Item index, relative to its position in the dataset.
 */
typedef int Index;

/** A cluster is simply a collection of indexes.
 *
 */
typedef std::vector<Index> Cluster;

/** A clustering is simply a collection of Clusters
 *
 */
typedef std::vector<Cluster> Clustering; 

////////////////////////////////////////////////////////////////////
struct Partition {
  
  /** A Label represents a cluster index, relative to its position in the clustering.
   */
  typedef int Label;
  
  /** An Index depicts an Item index, relative to its position in the dataset.
   */
  typedef std::map<Index,Label> Index2Label; 
  
  size_t nbrClusters() const { return m_labelSet.size(); }
  size_t nbrItems() const { return m_index2Label.size(); } 

  /// Converts itself to a clustering
  Clustering to_clustering() const {
    Clustering clustering( nbrClusters(), std::vector<int>() );
    for ( const auto& indexLabelPair: m_index2Label ) {  
    // for ( const std::pair<Index,Label>& indexLabelPair: m_index2Label ) {  
      clustering[ indexLabelPair.second ].push_back( indexLabelPair.first );
    }  
    return clustering;
  } 

  int cluster( int itemIdx ) const { return m_index2Label.at(itemIdx); } //@todo: getLabel +change Cluster -> Label
  void cluster( int itemIdx, int clusterIdx ) {  //@doto: setLabel
    m_index2Label[itemIdx] = clusterIdx; 
    m_labelSet.insert(clusterIdx);
  }

 private:
  std::set<int> m_labelSet; // A label is a cluster identifier. 
  Index2Label m_index2Label;
};

/** Convenient way to output the content of a clustering
 *
 */
inline std::ostream& operator<<( std::ostream& os, const Clustering& clt ) {
  os << "There are: " << clt.size() << " clusters\n";  
  for ( size_t c = 0; c < clt.size(); ++c ) {
    os << "cluster: " << c << std::endl;
    for ( size_t i = 0; i < clt.at(c).size() - 1; ++i ) {
      os << clt.at(c).at(i) << ", ";
    }
    os << clt.at(c).at(clt.at(c).size() - 1) << std::endl;
  }
  return os;  
}

/// Convenient way to output the content of a partition
inline std::ostream& operator<<( std::ostream& os, const Partition& p ) {
  return os << p.to_clustering();

}

/////////////////////////////////////////////////////////////////////////

} // namespace samogwas ends here. 

/****************************************************************************************/
#endif // SAMOGWAS_PARTITION_HPP 
