/****************************************************************************************
 * File: Clustering.hpp
 * Description: A structure representing a partition of the data set - i.e a division of the dataset into
 * -----------  non-overlapping subsets.
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 09/07/2014

 ***************************************************************************************/
#ifndef FLTM_CLUSTERING_HPP
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
  
  typedef std::map<int,int> Labels;
  
  size_t nbrClusters() const { return m_clusterSet.size(); }
  size_t nbrItems() const { return m_labels.size(); }

  /// Converts to a clustering and returns 
  Clustering to_clustering() const {
    Clustering clustering( nbrClusters(), std::vector<int>() );
    for ( const auto& i: m_labels ) {
      clustering[ i.second ].push_back( i.first );
    }  
    return clustering;
  } 

  int cluster( int item ) const { return m_labels.at(item); }
  void cluster( int item, int cluster ) {
    m_labels[item] = cluster;
    m_clusterSet.insert(cluster);
  }

 private:
  std::set<int> m_clusterSet;
  Labels m_labels;     
};

/** Convenient way to output the content of a clustering
 *
 */
inline std::ostream& operator<<( std::ostream& os, const Clustering& clt ) {
  os << "there are: " << clt.size() << " clusters\n";  
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

} // namespace fltmends here. fltm

/****************************************************************************************/
#endif // FLTM_CLUSTERING_HPP
