/****************************************************************************************
 * File: partition.hpp // CS I need to see the difference between Clustering.hpp and partition.hpp
 * Description: A data structure representing a partition of the data set - i.e a division of the data set into
 * -----------  non-overlapping subsets.
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
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

/** An Index depicts an Item index, relatively to its position in the dataset
 */
typedef int Index;

/** An Label represents a Label index, relatively to its position in whole clustering.
 */
typedef int Label;
/** A label can be seen as a cluster index
 *
 */

typedef int Label;

/** A cluster is simply a collection of Label (int)
 *
 */
typedef std::vector<Label> Cluster;

/** A clustering is simply a collection of Clusters
 *
 */
typedef std::vector<Cluster> Clustering; 

struct Partition {
  
  typedef std::map<Index,Label> Index2Label; // CS NOT VERY INFORMATIVE - should be named mapIntInt
  
  size_t nbrClusters() const { return m_clusterSet.size(); }
  size_t nbrItems() const { return m_index2Label.size(); } // CS nbrItems should be nbrOfKeys

  /// Converts itself to a clustering and returns // CS Converts WHAT into a clustering?
  Clustering to_clustering() const {
    Clustering clustering( nbrClusters(), std::vector<int>() );
    for ( const auto& i: m_index2Label ) { // CS INCOMPREHENSIBLE
      clustering[ i.second ].push_back( i.first );
    }  
    return clustering;
  } 

  int cluster( int item ) const { return m_index2Label.at(item); } // CS BAD IDENTIFIER, item should be named index (I presume)
  void cluster( int item, int cluster ) { // CS BAD IDENTIFIER clusterNumber
    m_index2Label[item] = cluster; // CS should be m_index2Label[index] = clusterNumber ???
    m_clusterSet.insert(cluster);
  }

 private:
  std::set<int> m_clusterSet;
  Index2Label m_index2Label; // CS     
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

/// Convenient way to output content of a partition
inline std::ostream& operator<<( std::ostream& os, const Partition& p ) {
  return os << p.to_clustering();

}

/////////////////////////////////////////////////////////////////////////

} // namespace samogwas ends here. S

/****************************************************************************************/
#endif // SAMOGWAS_PARTITION_HPP 
