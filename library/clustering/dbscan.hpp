/****************************************************************************************
 * File: DBSCAN.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 12/07/2014

 ***************************************************************************************/
#ifndef FLTM_DBSCAN_HPP
#define FLTM_DBSCAN_HPP

#include "clustering.hpp"
#include "partition.hpp"
// #include "distance/dissimilarity.hpp"

namespace samogwas
{

template<typename DistanceMatrix>
struct DBSCAN: public AlgoClust<DistanceMatrix> {
  typedef int Index;
  typedef std::vector<Index> Neighbors;
  typedef std::vector<Index> Labels;
  
  DBSCAN( DistanceMatrix* d, const int eles, const double eps ):
      AlgoClust<DistanceMatrix>(d), min_elems(eles), epsilon(eps) {   
  }
 
  virtual Partition operator()(); 
  virtual char* name() const {
    char* name = new char[80];
    sprintf( name, "DBSCAN_%d_%.3f", min_elems, epsilon);
    return name;
  }

 protected:
  Neighbors find_neighbors( const Index pid ) const;
  static Partition to_partition( const std::vector<int>& labels );
   
 protected:
  int min_elems;
  double epsilon;
};

} // namespace gwas ends here. fltm

/****************************************************************************************/
namespace samogwas
{

/**
 *
 */
template<typename DistanceMatrix>
Partition DBSCAN<DistanceMatrix>::operator()() { 
  size_t nvars = this->comp->size();
  std::vector<int> m_labels( nvars, -1);
  std::vector<int> visited( nvars, 0 );  
  int cluster_id = 0;
  for (int pid = 0; pid < nvars; ++pid) {
    if ( !visited[pid] ) {
      visited[pid] = 1;
      Neighbors ne = find_neighbors(pid);
      if ( ne.size() >= min_elems ) {
        m_labels[pid] = cluster_id; // partition.cluster( pid, cluster_id );
        for ( int i = 0; i < ne.size(); ++i) {
          int nPid = ne[i];
          if ( !visited[nPid] ) {
            visited[nPid] = 1;
            Neighbors ne1 = find_neighbors(nPid);
            if ( ne1.size() >= min_elems ) {
              for (const auto & n1 : ne1) {
                ne.push_back(n1);
              }
            }
          }
          if ( m_labels[nPid] == -1 ) {
            m_labels[nPid] = cluster_id;
          }
        }
        ++cluster_id;
      }
    }
  }

  return to_partition(m_labels);   
}

/**
 *
 */
template<typename DistanceMatrix>
typename DBSCAN<DistanceMatrix>::Neighbors DBSCAN<DistanceMatrix>::find_neighbors( const Index pid ) const {
  Neighbors ne;
  size_t nvars = this->comp->size();
  for ( Index i = 0; i < nvars; ++i ) {
    if ( this->comp->compute( i, pid ) <= epsilon ) {
      ne.push_back(i);
    }
  }

  return ne;
}

/////////////////////////////////////////
/**
 *
 */
template<typename DistanceMatrix>
Partition DBSCAN<DistanceMatrix>::to_partition( const std::vector<int>& labels ) {
  Partition partition;
  std::set<int> unique_labs;
  std::vector<int> singletons;
  for ( size_t i = 0; i < labels.size(); ++i ) {
    if ( labels.at(i) != -1 ) {
      partition.cluster( i,labels.at(i) );
      unique_labs.insert( labels.at(i) );
    } else {
      singletons.push_back(i);
    }
  }

  for ( auto& i: singletons ) {
    size_t cur_cluster = unique_labs.size();
    partition.cluster(i, cur_cluster);
    unique_labs.insert(cur_cluster);
  }
  
  return partition;
}



} // samogwas ends here

#endif // FLTM_DBSCAN_HPP
