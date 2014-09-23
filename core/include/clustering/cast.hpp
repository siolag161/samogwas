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
 * @ref: Ben Dor et al. 1999, Ben-Dor, Amir, Ron Shamir, and Zohar Yakhini (1999) Clustering gene expression patterns.
 * Journal of computational biology 6.3-4 : 281-297.
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 11/07/2014                                                       

 ***************************************************************************************/
#ifndef SAMOGWAS_CAST_HPP // CS Homogeneize the names all throughout the project.
#define SAMOGWAS_CAST_HPP

#include "distance/similarity.hpp"

#include "clustering.hpp" // AlgoClust
#include "partition.hpp"  // Partition, Clustering

namespace samogwas
{


/** The CAST methods returns a partition of the dataset, much like other types of AlgoClust.
 * Under the hood, during the computation, the private class CAST_Item is used to represent an 
 * item in the cluster. CAST_Item is not supposed to be used by an outside class.
 * CAST_Item provides a convenient way for holding, for each item, its current affinity relative to 
 * the current candidate cluster, and its (global) index relative to the original dataset.
 */
struct CAST_Item {  
  int globalIndex; 
  double affinity;
  CAST_Item( const int matIdx, const double aff): globalIndex(matIdx), affinity(aff) {}
};

/** The main class (functor). It takes as input a similarity matrix and a 
 *  threshold parameter thres (the same as t described above) in the constructor.
 */
template<typename SimiMatrix>
struct CAST: public AlgoClust<SimiMatrix>  {
  typedef std::vector<CAST_Item> CAST_Cluster;
  
  /** Constructor takes a Similarity matrix and a threshold criteron for removing and adding an item from the current cluster.
   *
   */
  CAST ( SimiMatrix* simMat, const double& thres): AlgoClust<SimiMatrix>(simMat), thresCAST(thres) {}

  /** Executes the algorithm given the current input.
   *
   */
  virtual Partition run() {
    /// Initializes the cluster grouping all the unassigned items with null affinity
    CAST_Cluster unAssignedCluster;
    for ( int i = 0; i < this->compMatrix->size(); ++i ) { // CS Nobody can guess that comp is a matrix. compMat should be used.
      unAssignedCluster.push_back(CAST_Item(i, 0.0) );    
    }
    /// Delegates the job to the method taking the newly initialized unassigned cluster.
    return run( unAssignedCluster );
  }

  /** Provides the name of the algorithm which its parameters
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
  
  // CS I do not see that there is such a step in the original CAST algorithm. @pdt: indeed, everytime a new openCluster is considered, gotta reset all the affinity relative to this one.
  // Instead I see : When a new cluster Copen is started, the initial affinities of all genes are 0 since Copen is empty.
  // CS Is it this step you describe here?
  // http://faculty.washington.edu/kayee/cluster/algs.pdf
  // Is it resetAffinity for the clusterofUnassignedObjects?
  
  /** Resets the current affinity matrix */
  void resetAffinity( CAST_Cluster& unAssignedCluster );

  
  // CS is unassignedCluster the equivalent of COpen? // pdt: not the same. COpen is the new open cluster to be formed, unassigned is a group of remaining objects
  // In CAST, objects are assigned to clusters.
  // It is incomprehensible that a cluster would be unassigned. Do you mean "under construction"?
  // Or is it clusterOfUnassignedObjects rather (as I guess)?
  void addGoodItem( CAST_Cluster& unassignedCluster, CAST_Cluster& openCluster, 
                    CompMatrix& simMatrix, const int itemIdx );

  /** Removes bad item from openCluster and moves them back to the unassignedCluster  */
  void removeBadItem( CAST_Cluster& unassignedCluster, CAST_Cluster& openCluster, 
                      CompMatrix& simMatrix, const int itemIdx );
  
  /** Moves the item indexed by itemIdx from sourceCluster, puts into targetCluster and then updates the affinity. */
  void updateClustersAffinity( CAST_Cluster& sourceCluster, CAST_Cluster& targetCluster, 
                               CompMatrix& simMatrix, const int itemIdx );
  
  /** Moves item indexed by itemIdx from sourceCluster to targetCluster  */
  void moveItemBetweenClusters( CAST_Cluster& source, CAST_Cluster& target, const int itemIdx );

 protected:
  // main parameter: CS unuseful 
  // CS the threshold involved in the computation of the CAST criterion for assigning an object to a cluster
  double thresCAST;     
};

//////////////////////////////////////////////////////////////////////////////
// CS What about the separation rules? Where to insert /////////// or /*********/?

/** Functor to find the index of the item of extremal affinity (mimimal or maximal).
 *
 */
struct ExtremumIndexAffinity { // CS heterogeneous notation: resetAffinity versus ExtremumIndexAffinity @pdt: not the same, one is a class, the other is a method
  template<typename Compare>
  const int operator()(const std::vector<CAST_Item>& clust, Compare comp) const;    
};

// CS What about the separation rules? Where to insert /////////// or /*********/?
}

/******************************************  IMPLEMENTATION BELOW **********************************************/
namespace samogwas
{


/** Implements the heuristic CAST algorithm following the original paper mentioned in the filehead.
 *  It requires a similarity matrix ( or any method that returns a similarity between a pair of item index)
 *  as input. It also requires a parameter t as the affinity threshold which determines whether we removes
 *  or add an item to the current candidate cluster. This parameter t also determines when a cluster is stabilizes.
 *  It returns a partition of the dataset
 */ 
template<typename SimiMatrix>
Partition CAST<SimiMatrix>::run( CAST_Cluster& unAssignedCluster ) {
  Partition result; // CS growingPartition ?
  while ( unAssignedCluster.size() ) {  // as long as there is still remains an object to be classified
    CAST_Cluster openCluster; // creates a new cluster
    resetAffinity( unAssignedCluster ); // reset the affinities of the objects remaining to be classified
    bool changesOccurred = true; 
    while (changesOccurred && unAssignedCluster.size()) { // begin organizing objets CS vague 
      // CS starts assigning objects to openCluster
      changesOccurred = false;
      ExtremumIndexAffinity maxCompute;
      ExtremumIndexAffinity minCompute;
      while ( unAssignedCluster.size() ) { // while there is still object to be assigned,
        // successively tries to assign the objects with the best affinities to openCluster
        int maxAffIdx = maxCompute( unAssignedCluster, std::greater<double>() );
        if ( unAssignedCluster.at(maxAffIdx).affinity >= thresCAST * openCluster.size() ) {
          changesOccurred = true;
          addGoodItem( unAssignedCluster, openCluster, *this->compMatrix, maxAffIdx );
        } else {
          break;
        }
      }
      while( unAssignedCluster.size() ) { // while there is still object to be assigned,
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
    int cluster_id =  result.nbrClusters(); // CS growingPartition ?
    for ( auto& item: openCluster ) { // CS INCOMPREHENSIBLE
      result.cluster( item.globalIndex, cluster_id ); // CS The procedure's name is ambiguous.
    }
  }  
  return result;  
}

//////////////////////////////////////////////////////////////////////////////////////
// sets all the affinity values to zeroes
template<typename SimiMatrix>
void CAST<SimiMatrix>::resetAffinity( CAST_Cluster& cluster ) { // CS resetAffinityValuesToZero
  for ( CAST_Cluster::iterator it = cluster.begin(); it != cluster.end(); ++it ) {
    it->affinity = 0.0;
  } 
}

////////////////////////////////////////////////////////////////////////////////////////

/// Delegates the task to updateClusterAffinity which moves the item indexed by itemIdx from unassignedCluster to openCluster and update the affinity correspondingly 
template<typename SimiMatrix>
void CAST<SimiMatrix>::addGoodItem( CAST_Cluster& unAssignedCluster, CAST_Cluster& openCluster, 
                                    CompMatrix& simiCompute, const int itemIdx ) {  
  updateClustersAffinity( unAssignedCluster, openCluster, 
                          simiCompute, itemIdx );  
}

/// Delegates the task to updateClusterAffinity which moves back the item indexed by itemIdx from openCluster to unassignedCluster and update the affinity correspondingly 
template<typename SimiMatrix>
void CAST<SimiMatrix>::removeBadItem( CAST_Cluster& unAssignedCluster, CAST_Cluster& openCluster, 
                                      CompMatrix& simiCompute, const int itemIdx ) {
  updateClustersAffinity( openCluster, unAssignedCluster, 
                          simiCompute, itemIdx); 
}

// Updates the affinity after changes are made to reflect current assignment
// CS updates the affinity values in sourceCluster and targetCluster after the move of indeXItem
template<typename SimiMatrix>
void CAST<SimiMatrix>::updateClustersAffinity( CAST_Cluster& sourceCluster, CAST_Cluster& targetCluster,
                                               CompMatrix& simiCompute, 
                                               const int itemIdx ) // CS bad denomination
    // indexOfItemInSourceCluster
{
  const int itemGlobalIndex = sourceCluster.at(itemIdx).globalIndex;
  moveItemBetweenClusters( sourceCluster, targetCluster, itemIdx ); 
  for (int i = 0; i < sourceCluster.size(); i++) {
    sourceCluster[i].affinity += simiCompute( sourceCluster[i].globalIndex, itemGlobalIndex ); // CS I am lost
    // Is type CompMatrix for ComparisonMatrix
    // or is it for simComputationFromMatrix?
    // suggestion for type identifier: CompSimFromMatrix? // pdt: it's a functor, used like a function
    // + itemIndex instead of item.globalIndex //pdt item.globalIndex is the index attribute from the item? 
  }

  for (int i = 0; i < targetCluster.size(); i++) {
    targetCluster[i].affinity += simiCompute( targetCluster[i].globalIndex, itemGlobalIndex ); // CS itemIndex 
  } 

}

/// Moves item from source cluster to target cluster
template<typename SimiMatrix>
    void CAST<SimiMatrix>::moveItemBetweenClusters( CAST_Cluster& source,
                                                    CAST_Cluster& target,
                                                    const int itemIdx ) { // indexOfItemInSourceCluster
      target.push_back( source.at(itemIdx) );
      //source.remove(itemIdx);
      source.erase( source.begin() + itemIdx ); // CS I do not understand this argument
      // CS C++ vector erase Removes the element at pos.
      // CS I would have expected indexOfItemInSourceCluster
    }

////////////////////////////////////////////////////////////////////////////////////////
// CS What is the rule for using //////////// or /*****************/?
/** Searches for the object in given cluster which has the optimal (best/worst) affinity
 *
 */
template<typename Compare>
const int ExtremumIndexAffinity::operator()( const std::vector<CAST_Item>& clust, Compare comp) const
    // CS obligatorily specify precondition
{ // CS change in presentation - left bracket should be after const
  int result = 0; // CS bad identifier idx_optimum
  for (int idx = 1; idx < clust.size(); idx++) {
    if ( comp(clust.at(idx).affinity, clust.at(result).affinity) ) { // CS not optimized
      result = idx;    
    }
  }
  
  return result;
}


}

#endif // SAMOGWAS_CAST_HPP
