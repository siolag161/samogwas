/****************************************************************************************
 * File: DBSCAN.hpp // CS incorrect
 * Description: An implementation of DBSCAN algorithm for clustering of data.h
 * // CS reference
 * 
 * The DBSCAN principle lies in constructing clusters from the estimated density distribution of the points to 
 * be clustered. This method requires two parameters: epsilon, the maximum radius of the neighborhood to be considered, 
 * and minPts, the minimum number of neighbor points needed for a cluster. DBSCAN exploits the fact that an point in a 
 * cluster also has its neighborhood in this cluster. The neighborhood of a point p is merely the set of 
 * points whose distance from p is less than or equal to epsilon. 
 * If this neighborhood is dense (i.e. its size is greater than or equal to minPts), then, p's neighborhood 
 * is grown as follows: 
 * through the addition of the proper neighborhoods of p'neighbors, provided that these neighborhoods are 
 * themselves dense (see lines 888 and 888). 
 * 
 * the neighborhoods of the non-visited neighbors of p are added to p's neighborhood, provided that these neighborhoods are 
 * themselves dense.
 * 
 * 
 * In the end, the grown neighborhood augmented with $o$ represents a new 
 * cluster. If the R-neighborhood of a point is not sufficiently dense, this point is labeled as noise (see line 888).
 * This point might later be found in the sufficiently dense R-neighborhood of a further visited point, and hence be 
 * assigned to the cluster constructed from this latter point (see line 888). Then, a new unvisited point is retrieved 
 * and processed, leading to the discovery of a further cluster or noise.\newline
 * 
 * @ref: Ester et al. (1996) Ester, M. and Kriegel, H.-P. and Sander, J. and Xu, X. A density-based algorithm for 
 * discovering clusters in large spatial databases with noise, Second International Conference on Knowledge Discovery 
 * and Data Mining (KDD96), 226-231.
 * 
 * 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet 
 * @date: 12/07/2014
 * 

 ***************************************************************************************/
#ifndef SAMOGWAS_DBSCAN_HPP 
#define SAMOGWAS_DBSCAN_HPP

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
                        // I do not see the usefulness of not keeping size_t as in other parts of the application. @pdt: it does not suppress size_t in other part. simply mean in this class we can use Index instead of size_t
  typedef std::vector<Index> Neighbors; // An important concept of DBSCAN involves Neighbors, a set of nearby points.
  typedef std::vector<Index> Labels; // A Labels maps an Item Index to its Cluster Index

  /** The constructor which takes three parameters: 
   * a distance matrix that holds the dissimilarities between pairs of point, 
   * the minimum number of a points to form a dense region, 
   * and the density threshold which determines whether two points are `close'.
   *  
   * 
   */
  DBSCAN( DistanceMatrix* d, const int eles, const double eps ):
      AlgoClust<DistanceMatrix>(d), minPts(eles), epsilon(eps) { // CS eles: not a nice identifier minElems (?)
                                                                    // CS eps: epsi would be more informative
  }

  /** The main method that executes the algorithm
   **/
  virtual Partition run();

  /// The method that returns the name of this algorithm, along with its parameters
  virtual char* name() const {
    char* name = new char[81];
    sprintf( name, "DBSCAN_%d_%.3f", minPts, epsilon);
    return name;
  }

 protected:
  /// Internal method that returns a set of `density-reachable' <points ???> from a given point.
  Neighbors find_neighbors( const Index pid ) const; // CS What is your rule for function identifiers all throughout 
                                                     // the application? findNeighbors versus find_neighbors

  /// Converts current clustering to the partition type
  static Partition toPartition( const Labels& labels ); // CS, in CAST, you use indexes, in DBSCAN, you use
                                                                   // labels? //@pdt: not the same thing
    
 protected:
  int minPts;
  double epsilon;
};

} // namespace gwas ends here. fltm CS ????

/****************************************************************************************/
namespace samogwas
{

/** The main method that executes the algorithm. The idea is based on the `density-reachability' concept. 
 * We visit every non-visited point in the data set and from there try to reach other non-visited points 
 * that are considered close to this one.
 * A point is reachable from another when the distance between them is less than the density-threshold parameter.
 * For this purpose, first we try to find for each point all of its direct neighbors. If the number of such neigbors
 * is below a given threshold minPts, we consider it as noisy point and move on to the next non-visited point.
 * Otherwise, we proceed to grow this group by trying to add to it other points that could be reach from other members of the group.
 * We continue until we cannot reach any more point outside of this group. A new cluster is then formed.
 * 
 */
template<typename DistanceMatrix> // CS bad type identifier, 
                                  // is it CompareMatrix as in clustering.hpp : AlgoClust( CompareMatrix* c): comp(c) {}
                                  // or is it really a dissimilarity matrix?
Partition DBSCAN<DistanceMatrix>::run() {
  size_t nvars = this->compMatrix->size(); // number of total variables
                                     // CS comp is not an informative identifier
  std::vector<size_t> m_labels( nvars, -1); // CS What is stored in m_labels?
  std::vector<size_t> visited( nvars, 0 ); // to keep track of visiting state for each point
  int cluster_id = 0; // initially there is no cluster formed  
  for (int pid = 0; pid < nvars; ++pid) { // CS we visit every non-visited point in this cluster
                                          // I do not see that there is a cluster, since you iterate on all points.
                                          // CS Why pid and not i?
                                          // CS use j or k
    if ( !visited[pid] ) {
      visited[pid] = 1;
      Neighbors neighbors  = find_neighbors(pid); // we find all the reachable points from the current point
                                          // heterogeneous style. finds all the reachable points from the current point
                                          // Vocabulary is very confusing: do not mix points and objets
                                          // Keep to one unique term.
                                          // + neighbs instead of ne
                                          
      if ( neighbors.size() >= minPts ) { // if the found region is not a noisy one, we form a new cluster
                                      // CS region???
                                      // CS + neighbs instead of ne
        m_labels[pid] = cluster_id; // partition.cluster( pid, cluster_id ); // we form a new cluster
                                    // CS bad identifier -> pointIdToClusterId ???
        for ( int i = 0; i < neighbors.size(); ++i) { // CS grows this cluster by trying to reach <CS other points?>
                                               // from each of its members
          int nPid = neighbors[i]; // 
          if ( !visited[nPid] ) { 
            visited[nPid] = 1;
            Neighbors subNeighbors = find_neighbors(nPid); // trying to find a new neighborhood // CS new neighbours
                                                  // CS neighbs1 
            if ( subNeighbors.size() >= minPts ) {
              for (const auto & neighbor : subNeighbors) { // CS Again, I hate this syntax. Is is necessary? Is it an optimization?
                neighbors.push_back(neighbor); // adds all the newly-found point to the current cluster
              }
            }
          }
          if ( m_labels[nPid] == -1 ) {
            m_labels[nPid] = cluster_id; // CS assigns point to its cluster <the cluster under construction>
          }
        }
        ++cluster_id; // increments the current cluster id
      }
    }
  }
  return toPartition(m_labels); // CS toPartition is a function?
}

/** A neighborhood of a given point is defined as all the points that are within a certain given distance.
 * CS The neighborhood of a given point is defined as all the points that are within a given distance.
 *
 */
template<typename DistanceMatrix>
typename DBSCAN<DistanceMatrix>::Neighbors DBSCAN<DistanceMatrix>::find_neighbors( const Index pid ) const {
  Neighbors ne;
  size_t nvars = this->compMatrix->size();
  for ( Index i = 0; i < nvars; ++i ) { // tries to add each point if the distance between it and the given point is within the threshold parameter 
    if ( this->compMatrix->compute( i, pid ) <= epsilon ) {
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
Partition DBSCAN<DistanceMatrix>::toPartition( const Labels& labels ) {
  Partition partition;
  std::set<int> unique_labs;
  std::vector<int> singletons;
  for ( size_t i = 0; i < labels.size(); ++i ) {
    if ( labels.at(i) != -1 ) {
      partition.cluster( i,labels.at(i) ); // CS identifier could be more informative.
      unique_labs.insert( labels.at(i) );
    } else {
      singletons.push_back(i);
    }
  }

  for ( auto& i: singletons ) {
    size_t cur_cluster = unique_labs.size();
    partition.cluster(i, cur_cluster);
    unique_labs.insert(cur_cluster); // CS semantics?
  }
  
  return partition;
}



} // samogwas ends here.

#endif // SAMOGWAS_DBSCAN_HPP
