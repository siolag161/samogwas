/****************************************************************************************
 * File: partition.hpp // CS I need to see the difference between Clustering.hpp and partition.hpp
 * Description: A data structure representing a partition of the data set - i.e a division of the data set into
 * -----------  non-overlapping subsets.
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 09/07/2014

 ***************************************************************************************/
#ifndef FLTM_CLUSTERING_HPP // CS (it should be SAMOGWAS_PARTITION_HPP)
#define FLTM_CLUSTERING_HPP

#include "distance/comparable.hpp"
#include <set>
#include <vector>
#include <map>

#include <iostream>
namespace samogwas
{

typedef std::vector< std::vector<int> > Clustering;

struct Partition {
  
  typedef std::map<int,int> Labels; // CS NOT VERY INFORMATIVE - should be named mapIntInt
  
  size_t nbrClusters() const { return m_clusterSet.size(); }
  size_t nbrItems() const { return m_labels.size(); } // CS nbrItems should be nbrOfKeys

  /// Converts to a clustering and returns // CS Converts WHAT into a clustering?
  Clustering to_clustering() const {
    Clustering clustering( nbrClusters(), std::vector<int>() );
    for ( const auto& i: m_labels ) { // CS INCOMPREHENSIBLE
      clustering[ i.second ].push_back( i.first );
    }  
    return clustering;
  } 

  int cluster( int item ) const { return m_labels.at(item); } // CS BAD IDENTIFIER, item should be named index (I presume)
  void cluster( int item, int cluster ) { // CS BAD IDENTIFIER clusterNumber
    m_labels[item] = cluster; // CS should be m_labels[index] = clusterNumber ???
    m_clusterSet.insert(cluster);
  }

 private:
  std::set<int> m_clusterSet;
  Labels m_labels; // CS     
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

} // namespace fltmends here. fltm // CS

/****************************************************************************************/
#endif // FLTM_CLUSTERING_HPP
