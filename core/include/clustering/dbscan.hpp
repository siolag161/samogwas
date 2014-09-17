/****************************************************************************************
 * File: DBSCAN.hpp // CS incorrect
 * Description: An implementation of DBSCAN algorithm for clustering of data.
 * // CS reference
 * @author: siolag161 (thanh.phan@outlook.com) // CS affiliations
 * @date: 12/07/2014

 ***************************************************************************************/
#ifndef FLTM_DBSCAN_HPP
#define FLTM_DBSCAN_HPP

#include "clustering.hpp"
#include "partition.hpp"

namespace samogwas
{

/** DBSCAN is a functor that implements AlgoClust interface for data clustering.
 * // CS Why not same commentary for CAST?
 *
 */
template<typename DistanceMatrix> 
// CS I refer to in clustering.hpp : AlgoClust( CompareMatrix* c): comp(c) {}
// CS For clarity, keep  a (good) type identifier throughout the whole application.
struct DBSCAN: public AlgoClust<DistanceMatrix> {
  
  typedef size_t Index; // CS I do not like Index as a type name
                        // I do not see the usefulness of not keeping size_t as in other parts of the application.    
  typedef std::vector<Index> Neighbors; // An important concept of DBSCAN involves Neighbors, a set of nearby objects.
  typedef std::vector<Index> Labels; 

  /** The constructor which takes three parameters: 
   * a distance matrix that holds the dissimilarities between pairs of object, 
   * the minimum number of a points to form a dense region, 
   * and the density threshold which determines whether two points are `close'.
   *  
   * 
   */
  DBSCAN( DistanceMatrix* d, const int eles, const double eps ):
      AlgoClust<DistanceMatrix>(d), min_elems(eles), epsilon(eps) { // CS eles: not a nice identifier minElems (?)
                                                                    // CS eps: epsi would be more informative
  }

  /** The main method that executes the algorithm
   **/
  virtual Partition run();

  /// The method that returns the name of this algorithm, along with its parameters
  virtual char* name() const {
    char* name = new char[81];
    sprintf( name, "DBSCAN_%d_%.3f", min_elems, epsilon);
    return name;
  }

 protected:
  /// Internal method that returns a set of `density-reachable' <points ???> from a given point.
  Neighbors find_neighbors( const Index pid ) const; // CS What is your rule for function identifiers all throughout 
                                                     // the application? findNeighbors versus find_neighbors

  /// Converts current clustering to the partition type
  static Partition to_partition( const std::vector<int>& labels ); // CS, in CAST, you use indexes, in DBSCAN, you use
                                                                   // labels?
                                                                   // to_partition versus toPartition?
    
 protected:
  int min_elems;
  double epsilon;
};

} // namespace gwas ends here. fltm CS ????

/****************************************************************************************/
namespace samogwas
{

/** The main method that executes the algorithm. The idea is based on the `density-reachability' concept. 
 * We visit every non-visited object in the data set and from there try to reach other non-visited objects 
 * that are considered close to this one. 
 * An object is reachable from another when the distance between them is less than the density-threshold parameter. 
 * A region <CS a region ??? a cluster ????> is formed once we cannot reach any more point outside of this group. If the 
 * cardinality of this region <CS group ?> is below the minPts <CS unknown> parameter, we discard <CS the cluster> and 
 * consider it as noise. Otherwise a new cluster is formed and we continue to proceed with the rest of non-visited objects ( if any b).
 *
 */
template<typename DistanceMatrix>
Partition DBSCAN<DistanceMatrix>::run() {
  size_t nvars = this->comp->size(); // number of total variables
  std::vector<int> m_labels( nvars, -1); 
  std::vector<int> visited( nvars, 0 ); // to keep track of visit state for every objects
  int cluster_id = 0; // initially there is no cluster formed
  printf("----------------------- DBSCAN: %d vars--------------------------\n", (int)nvars);
  
  for (int pid = 0; pid < nvars; ++pid) { // we visit every non-visited object in this cluster
    if ( !visited[pid] ) {
      visited[pid] = 1;
      Neighbors ne = find_neighbors(pid); // we find all the reachable objects from the current point
      if ( ne.size() >= min_elems ) { // if the found region is not a noisy one
        m_labels[pid] = cluster_id; // partition.cluster( pid, cluster_id ); // we form a new cluster
        for ( int i = 0; i < ne.size(); ++i) { // and proceed to grow this cluster by trying to reach from its memebers
          int nPid = ne[i]; // like above, we visit 
          if ( !visited[nPid] ) {
            visited[nPid] = 1;
            Neighbors ne1 = find_neighbors(nPid); // trying to find a new neighbors
            if ( ne1.size() >= min_elems ) {
              for (const auto & n1 : ne1) {
                ne.push_back(n1); // adds all the newly-found object to the current reachability region
              }
            }
          }
          if ( m_labels[nPid] == -1 ) {
            m_labels[nPid] = cluster_id; // assigns object to its cluster
          }
        }
        ++cluster_id; // increments the current cluster id
      }
    }
  }
  return to_partition(m_labels); //
}

/** A neighborhood of a given object is defined as all the points that are within a certain given distance.
 *
 */
template<typename DistanceMatrix>
typename DBSCAN<DistanceMatrix>::Neighbors DBSCAN<DistanceMatrix>::find_neighbors( const Index pid ) const {
  Neighbors ne;
  size_t nvars = this->comp->size();
  for ( Index i = 0; i < nvars; ++i ) { // Simply try to add every point if the distance between it and the given object is within the threshold parameter 
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
