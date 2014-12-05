/****************************************************************************************
 * File: cast.hpp
 * Description: CAST (Clustering Affinity Search Technique) is a clustering algorithm proposed in Ben Dor et al. 1999 
 * The idea of the algorithm is to represent the true underlying data as a probabilistic graph. A cluster 
 * is modeled as a clique of this graph. 
 * 
 *The implemention below is a heuristic version described in the original paper. The idea is while we still have 
 * items that remain to be clustered (grouped together in a set called unassignedCluster), we successively try 
 * to generate a new cluster (openCluster). 
 * 
 * The affinity measure of an object is defined as the average of the similarity measures between objects in 
 * a given cluster and another object outside this cluster. To grow a cluster, the CAST algorithm successively 
 * adds the object with greatest affinity to this cluster as long as this maximum affinity satisfies a threshold 
 * constraint (t). When the maximum affinity drops below the threshold, CAST removes the object with the minimum affinity 
 * with respect to the cluster. Additions and removals are operated until the new cluster (openCluster) stabilizes. 
 * Finally, the cluster is closed. 
 *
 * @ref: Ben Dor et al. (1999) Ben-Dor, A., Ron S., and Zohar Y. Clustering gene expression patterns.
 * Journal of computational biology, 6(3-4), 281-297.
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 11/07/2014                                                       

 ***************************************************************************************/
#ifndef SAMOGWAS_CAST_HPP // CS Homogeneize the names all throughout the project.
#define SAMOGWAS_CAST_HPP

#include <memory>

#include "distance/similarity.hpp" // SimilarityMatrix

#include "clustering.hpp" // AlgoClust
#include "partition.hpp"  // Partition, Clustering

namespace samogwas
{


/** The CAST method returns a partition of the dataset, much like other types of AlgoClust.
 * Under the hood, during the computation, the private class CAST_Item is used to represent an 
 * item in the cluster. CAST_Item is not supposed to be used by an outside class.
 * CAST_Item provides a convenient way for holding, for each item, its current affinity relative to 
 * the current candidate Opencluster, and its (global) index relative to the original dataset.
 */
struct CAST_Item {  
  int globalIndex; 
  double affinity;
  CAST_Item( const int matIdx, const double aff): globalIndex(matIdx), affinity(aff) {}
};

/** CAST is a functor that implements AlgoClust interface for data clustering.
 *  It takes as input a similarity matrix and a 
 *  threshold parameter thres (the same as t described above) in the constructor.
 */
template<typename SimiMatrix> //  for instance, may be instanciated as MutInfoSimilarityMatrix
struct CAST: public AlgoClust<SimiMatrix>  {
  typedef std::vector<CAST_Item> CAST_Cluster;
  
  /** Constructor takes a similarity matrix and a threshold criterion for removing and adding an item from the current 
   * candidate Opencluster.
   *
   */
  CAST ( std::shared_ptr<SimiMatrix> simMat, const double& thres): AlgoClust<SimiMatrix>(simMat), thresCAST(thres) {}

  /** Executes the algorithm given the current input (SimiMatrix).
   *
   */
  virtual Partition run() {
    /// Initializes the cluster unAssignedCluster, grouping all the unassigned items.
    CAST_Cluster unAssignedCluster; 
    for ( int i = 0; i < this->compMatrix->nbrVariables(); ++i ) { // SimilarityMatrix has type SimiMatrix (for instance
                                                           // MutInfoSimilarityMatrix)
                                                           // @todo: remove direct access, and change the name
      unAssignedCluster.push_back(CAST_Item(i, 0.0) );    
    }
    /// Delegates the job to the method taking the newly initialized unassignedCluster as parameter.
    return run( unAssignedCluster );
  }

  /** Provides the name of the algorithm along with its parameters
   *
   */
  virtual char* name() const {
    char* name = new char[80];
    sprintf( name, "CAST_%.3f", thresCAST);
    return name;
  }
  
 protected:
  /** Performs on an unAssignedCluster, setup by the run() method above
   */
  inline Partition run( CAST_Cluster& unassignedCluster );
  
  
  /** Resets the affinities of the items in unAssignedCluster */
  void resetAffinities( CAST_Cluster& unAssignedCluster );

  
  // Adds good item to OpenCluster and removes it from unassignedCluster 
  void addGoodItem( CAST_Cluster& unassignedCluster, CAST_Cluster& openCluster, 
                    SimilarityMatrix& simMatrix, const int itemIdx );

  /** Removes bad item from openCluster and moves it back to unassignedCluster  */
  void removeBadItem( CAST_Cluster& unassignedCluster, CAST_Cluster& openCluster, 
                      SimilarityMatrix& simMatrix, const int itemIdx );
  
  /** Moves the item indexed by itemIdx from sourceCluster, puts into targetCluster and then updates the affinity. */
  void updateClustersAffinity( CAST_Cluster& sourceCluster, CAST_Cluster& targetCluster, 
                               SimilarityMatrix& simMatrix, const int itemIdx );
  
  /** Moves item indexed by itemIdx from sourceCluster to targetCluster  */
  void moveItemBetweenClusters( CAST_Cluster& source, CAST_Cluster& target, const int itemIdx );

 protected:
  // main parameter of CAST
  double thresCAST;     
};

//////////////////////////////////////////////////////////////////////////////

/** Functor to find the index of the item of extremal affinity (mimimal or maximal).
 *
 */
struct ExtremumIndexAffinity {
  template<typename Compare>
  const int operator()(const std::vector<CAST_Item>& cluster, Compare comparat) const;    
};

}

/******************************************  IMPLEMENTATION BELOW **********************************************/
namespace samogwas
{


/** Implements the heuristic CAST algorithm following the original paper mentioned in the header.
 *  This heuristic requires as input a similarity matrix ( or any method that returns a similarity between 
 *  a pair of item indexes).
 .  It also requires a parameter t as the affinity threshold which determines whether we remove
 *  or add an item to the current candidate Opencluster. This parameter t also determines when a cluster is stabilized.
 *  This method returns a partition of the dataset.
 */ 
template<typename SimiMatrix> 
Partition CAST<SimiMatrix>::run( CAST_Cluster& unAssignedCluster ) {
  Partition result; 
  while ( unAssignedCluster.size() ) {  // as long as there still remains an object to be classified (i.e. unassigned)
    CAST_Cluster openCluster; // creates a new cluster
    resetAffinities( unAssignedCluster ); // reset the affinities of the unassigned objects
    bool changesOccurred = true; 
    while (changesOccurred && unAssignedCluster.size()) { // starts assigning objects to openCluster
      changesOccurred = false;
      ExtremumIndexAffinity maxCompute;
      ExtremumIndexAffinity minCompute;
      while ( unAssignedCluster.size() ) { // while there is still an object to be assigned,
        // successively tries to assign the objects with the best affinities to openCluster
        int maxAffIdx = maxCompute( unAssignedCluster, std::greater<double>() );
        if ( unAssignedCluster.at(maxAffIdx).affinity >= thresCAST * openCluster.size() ) {
          changesOccurred = true;
          addGoodItem( unAssignedCluster, openCluster, *this->compMatrix, maxAffIdx );
        } else {
          break;
        }
      }
      while( unAssignedCluster.size() ) { // while there is still an object to be assigned,
        // successively tries to remove the objects with the worst affinities from openCluster
        int minAffIdx = minCompute( openCluster, std::less<double>() );
        if ( openCluster.at(minAffIdx).affinity < thresCAST * openCluster.size() ) {
          changesOccurred = true;
          removeBadItem( unAssignedCluster, openCluster, *this->compMatrix, minAffIdx );
        } else {
          break;
        }       
      }
    }  // until stabilized
    // puts the newly created openCluster into the current partition and continues
    int cluster_id =  result.nbrClusters(); 
    for ( auto& item: openCluster ) {
      result.setLabel( item.globalIndex, cluster_id );
    }
  }  
  return result;  
}

//////////////////////////////////////////////////////////////////////////////////////
// Sets all the affinity values to zeroes.
template<typename SimiMatrix>
void CAST<SimiMatrix>::resetAffinities( CAST_Cluster& cluster ) {
  for ( CAST_Cluster::iterator it = cluster.begin(); it != cluster.end(); ++it ) {
    it->affinity = 0.0;
  } 
}

//////////////////////////////////////////////////////////////////////////////////////
// Delegates the task to updateClusterAffinity which moves the item indexed by itemIdx from unassignedCluster 
// to openCluster and update the affinity correspondingly. 
template<typename SimiMatrix>
void CAST<SimiMatrix>::addGoodItem( CAST_Cluster& unAssignedCluster, CAST_Cluster& openCluster, 
                                    SimilarityMatrix& simiCompute, const int itemIdx ) {  
  updateClustersAffinity( unAssignedCluster, openCluster, 
                          simiCompute, itemIdx );  
}

//////////////////////////////////////////////////////////////////////////////////////
// Delegates the task to updateClusterAffinity which moves back the item indexed by itemIdx from openCluster to 
// unassignedCluster and update the affinity correspondingly 
template<typename SimiMatrix>
void CAST<SimiMatrix>::removeBadItem( CAST_Cluster& unAssignedCluster, CAST_Cluster& openCluster, 
                                      SimilarityMatrix& simiCompute, const int itemIdx ) {
  updateClustersAffinity( openCluster, unAssignedCluster, 
                          simiCompute, itemIdx); 
}

//////////////////////////////////////////////////////////////////////////////////////
// Updates the affinity values in sourceCluster and targetCluster after the move of indeXItem.
template<typename SimiMatrix>
void CAST<SimiMatrix>::updateClustersAffinity( CAST_Cluster& sourceCluster, CAST_Cluster& targetCluster,
                                               SimilarityMatrix& simiCompute, 
                                               const int itemIdx )
{
  const int itemGlobalIndex = sourceCluster.at(itemIdx).globalIndex;
  moveItemBetweenClusters( sourceCluster, targetCluster, itemIdx ); 
  for (int i = 0; i < sourceCluster.size(); i++) {
    sourceCluster[i].affinity += simiCompute( sourceCluster[i].globalIndex, itemGlobalIndex ); 
  }

  for (int i = 0; i < targetCluster.size(); i++) {
    targetCluster[i].affinity += simiCompute( targetCluster[i].globalIndex, itemGlobalIndex ); 
  } 

}

//////////////////////////////////////////////////////////////////////////////////////
// Moves item from source cluster to target cluster.
template<typename SimiMatrix>
void CAST<SimiMatrix>::moveItemBetweenClusters( CAST_Cluster& source,
                                                CAST_Cluster& target,
        const int itemIdx ) {
 target.push_back( source.at(itemIdx) );
 source.erase( source.begin() + itemIdx );
}

////////////////////////////////////////////////////////////////////////////////////////
/** Searches for the object in given cluster which has the optimal (best/worst) affinity.
 *
 */
template<typename Compare>
const int ExtremumIndexAffinity::operator()( const std::vector<CAST_Item>& cluster, Compare comparat) const
{ 
  int result = 0; 
  for (int idx = 1; idx < cluster.size(); idx++) {
    if ( comparat(cluster.at(idx).affinity, cluster.at(result).affinity) ) {
      result = idx;    
    }
  }
  
  return result;
}


} // namespace samogwas ends here.

#endif // SAMOGWAS_CAST_HPP
