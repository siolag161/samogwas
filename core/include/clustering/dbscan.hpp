/****************************************************************************************
 * File: dbscan.hpp
 * Description: An implementation of the DBSCAN algorithm dedicated to data clustering (Ester et al. (1996)).
 * 
 * 
 * The DBSCAN principle lies in constructing clusters from the estimated density distribution of the points to 
 * be clustered. This method requires two parameters: epsilon, the maximum radius of the neighborhood to be considered, 
 * and minPts, the minimum number of neighbor points needed for a cluster. DBSCAN exploits the fact that a point in a 
 * cluster also has its neighborhood in this cluster. The neighborhood of a point p is merely the set of 
 * points whose distance from p is less than or equal to epsilon. 
 * 
 * If this neighborhood is dense (i.e. its size is greater than or equal to minPts), then, p's neighborhood 
 * is grown as follows: 
 * the proper neighborhoods of p's neighbors are successively added to p's growing neighborhood. The condition to add a 
 * neighborhood is that it is itself dense. Thus, the newly added point is considered 
 * the neighbor of p. This growing process is iterated until there is no more point left unvisited in p's neighborhood. Thus, 
 * a new cluster has been created.
 * 
 * Otherwise (if p's neighborhood is not dense), the next unvisited point is considered. (The point p which is skipped
 * may further be integrated in a cluster.)
 * 
 * The global process is terminated when there is no more unvisited point left. 
 * The DBSCAN clustering algorithm is likely to output points unassigned to clusters (i.e. noisy points).
 * 
 * @ref: Ester et al. (1996) Ester, M. and Kriegel, H.-P. and Sander, J. and Xu, X. A density-based algorithm for 
 * discovering clusters in large spatial databases with noise, Second International Conference on Knowledge Discovery 
 * and Data Mining (KDD96), 226-231.
 * 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet 
 * @date: 12/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_DBSCAN_HPP 
#define SAMOGWAS_DBSCAN_HPP

#include <memory>

#include "clustering.hpp" // AlgoClust
#include "partition.hpp" // Partition

namespace samogwas
{

/** DBSCAN is a functor that implements AlgoClust interface for data clustering.
*   The DBSCAN method returns a partition of the dataset, much like other types of AlgoClust.
*/
 
template<typename DissMatrix>  // for instance, may be instanciated as MutInfoDissimilarityMatrix
struct DBSCAN: public AlgoClust<DissMatrix> {
  enum { UNASSIGNED_LABEL = -1 };
  typedef int Index;
  typedef std::vector<Index> Neighbors; // An important concept of DBSCAN involves Neighbors, a set of nearby points.typedef int Index;
    
  typedef int Label;
  typedef std::vector<Index> LabelSet; // A Label maps a point Index to its cluster Index.
                                     // A  for instance, may be instanciated as MutInfoSimilarityMatrix label indicates that the point is not yet assigned to any cluster.
  /** The constructor takes three parameters: 
   * a distance matrix that holds the dissimilarities between pairs of points, 
   * the minimum number of points to form a dense neighborhood, 
   * and the maximum radius (epsilon) of the neighborhood to be considered, which determines whether two points are neighbors.
   */
   
  DBSCAN( std::shared_ptr<DissMatrix> d, const int minPoints, const double epsi ):
      AlgoClust<DissMatrix>(d), minPts(minPoints), epsilon(epsi) { 
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
  /// Internal method that returns the set of neighbors for a given point.
  Neighbors find_neighbors( const Index pid ) const; 

  /// Converts current clustering to the partition type
  static Partition toPartition( const LabelSet& LabelSet ); 
    
 protected:
  int minPts;
  double epsilon;
};

} // namespace samogwas ends here.

/************************** IMPLEMENTATION BELOW **************************************************************/
namespace samogwas
{

/** The main method that executes the algorithm. The idea is based on the `density-reachability' concept. 
 * We visit every non-visited point in the data set and from there try to reach other non-visited points 
 * that are close to this one.
 * A point is reachable from another one when the dissimilarity between them is less than the radius 
 * parameter (epsilon).
 * For this purpose, for each non-visited point, we first identify its direct neighbors. If the number of such neigbors
 * is below a given threshold minPts, we consider it as noisy point and move on to the next non-visited point.
 * Otherwise, we proceed to grow this neighborhood by trying to add to it other points that could be reached from other 
 * members of the neighborhood.
 * We continue until we cannot reach any more point outside of this group. A new cluster is then formed.
 */
template<typename DissMatrix> 
Partition DBSCAN<DissMatrix>::run() {

  size_t nvars = this->compMatrix->nbrVariables(); // total number of variables
  LabelSet m_LabelSet( nvars,  UNASSIGNED_LABEL); 
  std::vector<size_t> visited( nvars, 0 ); // to keep track of visiting state for each point
  int cluster_id = 0; // Initially there is no cluster formed.
  for (int pid = 0; pid < nvars; ++pid) { // We visit every non-visited point. 
                                         
    if ( !visited[pid] ) {
      visited[pid] = 1;
      Neighbors neighbors  = find_neighbors(pid); 
      if ( neighbors.size() >= minPts ) { // If the neighborhood is dense,
        m_LabelSet[pid] = cluster_id; // we form a new cluster.
        // printf("clust(%d): %d\n", pid, cluster_id);
        for ( int i = 0; i < neighbors.size(); ++i) { // We grow this cluster by trying to reach other points
                                                     // from each of its members.
          int nPid = neighbors[i]; // 
          if ( !visited[nPid] ) { 
            visited[nPid] = 1;
            Neighbors subNeighbors = find_neighbors(nPid); // trying to find a new dense neighborhood
            if ( subNeighbors.size() >= minPts ) {
              for (const auto & neighbor : subNeighbors) { 
                neighbors.push_back(neighbor); // adds all the newly found points to the current cluster
              }
            }
          }
          if ( m_LabelSet[nPid] ==  UNASSIGNED_LABEL ) { // to avoid overriding a possible previous cluster assignment
                    // printf("clust(%d): %d\n", nPid, cluster_id);

            m_LabelSet[nPid] = cluster_id; 
          }
        }
        ++cluster_id; // increments the current cluster id
      }
    }
  }
  return toPartition(m_LabelSet); 
}

/////////////////////////////////////////
/** A neighborhood of a given point is defined as all the points that are within a certain given radius (epsilon).
 */
template<typename DissMatrix>
typename DBSCAN<DissMatrix>::Neighbors DBSCAN<DissMatrix>::find_neighbors( const Index pid ) const {
  Neighbors ne;
  size_t nvars = this->compMatrix->nbrVariables(); // @todo: remove direct access to compMatrix
  for ( Index i = 0; i < nvars; ++i ) {
    // printf("diff(%d,%d) = %f vs %f\n", pid, i, this->compMatrix->compute( i, pid ), epsilon);

    if ( this->compMatrix->compute( i, pid ) <= epsilon ) { // @todo: remove direct access to compMatrix
      
      ne.push_back(i);
    }

  }

  return ne;
}

/////////////////////////////////////////
template<typename DissMatrix>
Partition DBSCAN<DissMatrix>::toPartition( const LabelSet& LabelSet ) {
  Partition partition;
  std::set<Label> labs;
  std::vector<int> singletons;
  for ( size_t i = 0; i < LabelSet.size(); ++i ) {
    if ( LabelSet.at(i) != UNASSIGNED_LABEL ) {
      partition.setLabel( i, LabelSet.at(i) ); // CS identifier could be more informative.
                                           // The label is the cluster identifier.
      labs.insert( LabelSet.at(i) );
    } else {
      singletons.push_back(i);
    }
  }

  for ( auto& i: singletons ) {
    size_t curCluster = labs.size();
    partition.setLabel(i, curCluster);
    labs.insert(curCluster);
  }
  
  return partition;
}



} // samogwas ends here.

#endif // SAMOGWAS_DBSCAN_HPP
