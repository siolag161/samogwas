/****************************************************************************************
 * File: CAST.hpp
 * Description: CAST is an algorith for clustering which is based on the usage of a
 * -----------  dissimilarity matrix. 
 * // CS This description could apply to other algorithms.
 * // At least, cite the publication and technical report.
 * @author: siolag161 (thanh.phan@outlook.com) // CS Remark about affiliations not taken into account, since
 *                                             // the first checking round in autumn 2013 
 * @date: 11/07/2014

 ***************************************************************************************/
#ifndef CLUSTERING_CAST_HPP // CS Homogeneize the names all throughout the project.
#define CLUSTERING_CAST_HPP

#include "distance/similarity.hpp"

#include "clustering.hpp"
#include "partition.hpp"

namespace samogwas
{


/** The CAST methods returns a repartition CS repartition or partition ??? of the dataset 
 * CS CS CLARIFY in the form of non-overlapping set sets ??? of <objects of type> CAST_Item ???.
 *  A CAST_Item is a structure which holds an index of the object ( <CS suggestion index relative to the> relatively in the original dataset ) and
 * the current object's affinity computed relatively to its current assigned cluster.
 *
 */
struct CAST_Item {  
  int index;
  double affinity;
  CAST_Item(const int matIdx, // const int glbIdx, // CS Why this specification?
            const double aff):
      index(matIdx), // globalIndex(glbIdx)
      affinity(aff) {}
};
////////////////////////////////////////////////////////////////////////////
template<typename SimiMatrix>
struct CAST: public AlgoClust<SimiMatrix>  {
  /**
   *
   */
  CAST ( SimiMatrix* sim, const double& thres): AlgoClust<SimiMatrix>(sim), thresCAST(thres) {}
  // CS It is not the first time I demand that sim be renamed in the more informative name simMat.

  /**
   *
   */
  virtual Partition run() {
    std::vector<CAST_Item> unAssignedCluster;
    for ( int i = 0; i < this->comp->size(); ++i ) { // CS Nobody can guess that comp is a matrix. compMat should be used.
      unAssignedCluster.push_back(CAST_Item(i, 0.0) );    
    }
    return run( unAssignedCluster );
  }

  /** Provide the name of the algorithm which its parameters
   *
   */
  virtual char* name() const {
    char* name = new char[80];
    sprintf( name, "CAST_%.3f", thresCAST);
    return name;
  }
  
 protected:
  /** Performs on an unAssignedCluster, setup by the run() method above
   *
   */
  inline Partition run( std::vector<CAST_Item>& unassignedCluster );

  /** Resets the current affinity matrix
   *
   */
  inline void resetAffinity( std::vector<CAST_Item>& unAssignedCluster );

 protected:
  // main parameter: CS unuseful 
  // CS the threshold involved in the computation of the CAST criterion for assigning an object to a cluster
  double thresCAST;     
};

//////////////////////////////////////////////////////////////////////////////

struct AffinityCompute { // CS heterogeneous notation: resetAffinity versus AffinityCompute 
  template<typename Compare>
  const int operator()(const std::vector<CAST_Item>& clust, Compare comp) const;    
};

/** Resets the current affinity of the cluster given as parameter
 *
 */
inline void resetAffinity( std::vector<CAST_Item>& cluster);
// What is the link with the above inline void resetAffinity?
//
// CS I do not see that there is such a step in the original CAST algorithm.
// Instead I see : When a new cluster Copen is started, the initial affinities of all genes are 0 since Copen is empty.
// CS Is it this step you describe here?
// http://faculty.washington.edu/kayee/cluster/algs.pdf
// Is it resetAffinity for the clusterofUnassignedObjects?


// CS is unassignedCluster the equivalent of COpen?
// In CAST, objects are assigned to clusters.
// It is incomprehensible that a cluster would be unassigned. Do you mean "under construction"?
// Or is it clusterOfUnassignedObjects rather (as I guess)?
inline void addGoodItem( std::vector<CAST_Item>& unassignedCluster, std::vector<CAST_Item>& openCluster, 
                         CompMatrix& simMatrix, const int clusterIdx );

inline void removeBadItem( std::vector<CAST_Item>& unassignedCluster, std::vector<CAST_Item>& openCluster, 
                           CompMatrix& simMatrix, const int clusterIdx );

inline void updateClustersAffinity( std::vector<CAST_Item>& sourceCluster, std::vector<CAST_Item>& targetCluster, 
                                    CompMatrix& simMatrix,
                                    const int clusterIndex );
// CS heterogeneous notation: clusterIdx versus clusterIndex                                    

inline void moveItemBetweenClusters( std::vector<CAST_Item>& source, std::vector<CAST_Item>& target, const int clusterIndex );

}

/****************************************************************************************/
namespace samogwas
{


/** Implements the CAST algorithm following the original paper mentioned in the filehead
 *
 */
template<typename SimiMatrix>
Partition CAST<SimiMatrix>::run( std::vector<CAST_Item>& unAssignedCluster ) {
  Partition result;
  while ( unAssignedCluster.size() ) {  // as long as there is still remains an object to be classified
    std::vector<CAST_Item> openCluster; // creates a new cluster
    resetAffinity( unAssignedCluster ); // reset the affinities of the objects remaining to be classified
    bool changesOccurred = true; 
    while (changesOccurred && unAssignedCluster.size()) { // begin organizing objets CS vague 
                                                          // CS starts assigning objects to openCluster
      changesOccurred = false;
      AffinityCompute maxCompute;
      AffinityCompute minCompute;
      while ( unAssignedCluster.size() ) { // tries to put the best-affinity objects in the open cluster
                                           // CS tries to assign the objects with the best affinities to openCluster
        int maxAffIdx = maxCompute( unAssignedCluster, std::greater<double>() );
        if ( unAssignedCluster.at(maxAffIdx).affinity >= thresCAST*openCluster.size() ) {
          changesOccurred = true;
          addGoodItem( unAssignedCluster, openCluster, *this->comp, maxAffIdx );
        } else {
          break;
        }
      }
      while( unAssignedCluster.size() ) { // tries to remove the worst-affinity objects from the open cluster
        int minAffIdx = minCompute( openCluster, std::less<double>() );
        if ( openCluster.at(minAffIdx).affinity < thresCAST*openCluster.size() ) {
          changesOccurred = true;
          removeBadItem( unAssignedCluster, openCluster, *this->comp, minAffIdx );
        } else {
          break;
        }       
      }
    }  // until stablized
    // put the newly created open cluster into the current repartition and continues
    int cluster_id =  result.nbrClusters(); 
    for ( auto& item: openCluster ) {
      result.cluster( item.index, cluster_id );
    }
  }  
  return result;  
}

// sets all the affinity valuse to zeroes
template<typename SimiMatrix>
void CAST<SimiMatrix>::resetAffinity( std::vector<CAST_Item>& cluster ) {
  for ( std::vector<CAST_Item>::iterator it = cluster.begin(); it != cluster.end(); ++it ) {
    it->affinity = 0.0;
  } 
}

////////////////////////////////////////////////////////////////////////////////////////
/** Searches for the object in given cluster which has the optimal (best/worst) affinity
 *
 */
template<typename Compare>
const int AffinityCompute::operator()( const std::vector<CAST_Item>& clust, Compare comp) const
{
  int result = 0; 
  for (int idx = 1; idx < clust.size(); idx++) {
    if ( comp(clust.at(idx).affinity, clust.at(result).affinity) ) { 
      result = idx;    
    }
  }
  
  return result;
}

////////////////////////////////////////////////////////////////////////////////////////

/// Adds a candidate item to the open cluster
void addGoodItem( std::vector<CAST_Item>& unAssignedCluster, std::vector<CAST_Item>& openCluster, 
                  CompMatrix& simCompute, const int clusterIdx ) {  
  updateClustersAffinity( unAssignedCluster, openCluster, simCompute,
                          clusterIdx );  
}

/// removes a potentially bad item from the cluster
void removeBadItem( std::vector<CAST_Item>& unAssignedCluster, std::vector<CAST_Item>& openCluster, 
                    CompMatrix& simCompute, const int clusterIdx ) {
  updateClustersAffinity( openCluster, unAssignedCluster, simCompute,
                          clusterIdx); 
}

/// updates the affinity after changes are made to reflect current repartition
void updateClustersAffinity( std::vector<CAST_Item>& sourceCluster, std::vector<CAST_Item>& targetCluster,
                             CompMatrix& simCompute, 
                             const int clusterIndex )
{
  const CAST_Item item = sourceCluster.at(clusterIndex);
  moveItemBetweenClusters( sourceCluster, targetCluster, clusterIndex );
  for (int i = 0; i < sourceCluster.size(); i++) {
    sourceCluster[i].affinity += simCompute( sourceCluster[i].index, item.index );
  }

  for (int i = 0; i < targetCluster.size(); i++) {
    targetCluster[i].affinity += simCompute( targetCluster[i].index, item.index );
  } 

}

/// Moves item from source cluster to target cluster
void moveItemBetweenClusters( std::vector<CAST_Item>& source,
                              std::vector<CAST_Item>& target,
                              const int clusterIndex ) {
  target.push_back( source.at(clusterIndex) );
  //source.remove(clusterIndex);
  source.erase( source.begin() + clusterIndex );
}


}

#endif // CLUSTERING_CAST_HPP
