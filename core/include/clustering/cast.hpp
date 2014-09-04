/****************************************************************************************
 * File: CAST.hpp
 * Description: CAST is an algorith for clustering which is based on the usage of a
 * -----------  dissimilarity matrix. 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 11/07/2014

 ***************************************************************************************/
#ifndef CLUSTERING_CAST_HPP
#define CLUSTERING_CAST_HPP

#include "distance/similarity.hpp"

#include "clustering.hpp"
#include "partition.hpp"

namespace samogwas
{


/** The CAST methods returns a repartition of the dataset in the form of non-overlapping set of CAST_Item.
 *  A CAST_Item is a structure which holds an index of of the object ( relatively in the original dataset ) and
 * the current affinity computed relatively to its current cluster.
 *
 */
struct CAST_Item {  
  int index;
  double affinity;
  CAST_Item(const int matIdx, // const int glbIdx,
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

  /**
   *
   */
  virtual Partition run() {
    std::vector<CAST_Item> unAssignedCluster;
    for ( int i = 0; i < this->comp->size(); ++i ) { 
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

  /** Reset the current affinity matrix
   *
   */
  inline void resetAffinity( std::vector<CAST_Item>& unAssignedCluster );

 protected:
  // main parameter: the threshold for assigning objects
  double thresCAST;     
};

//////////////////////////////////////////////////////////////////////////////
struct AffinityCompute {
  template<typename Compare>
  const int operator()(const std::vector<CAST_Item>& clust, Compare comp) const;    
};

/** Resets the current affinity of the given cluster as parameter
 *
 */
inline void resetAffinity( std::vector<CAST_Item>& cluster);


inline void addGoodItem( std::vector<CAST_Item>& unassignedCluster, std::vector<CAST_Item>& openCluster, 
                         CompMatrix& simMatrix, const int clusterIdx );

inline void removeBadItem( std::vector<CAST_Item>& unassignedCluster, std::vector<CAST_Item>& openCluster, 
                           CompMatrix& simMatrix, const int clusterIdx );

inline void updateClustersAffinity( std::vector<CAST_Item>& sourceCluster, std::vector<CAST_Item>& targetCluster, 
                                    CompMatrix& simMatrix,
                                    const int clusterIndex );

inline void moveItemBetweenClusters( std::vector<CAST_Item>& source, std::vector<CAST_Item>& target, const int clusterIndex );

}

/****************************************************************************************/
namespace samogwas
{


/** Implements the CAST algorithm following the original paper
 *
 */
template<typename SimiMatrix>
Partition CAST<SimiMatrix>::run( std::vector<CAST_Item>& unAssignedCluster ) {
  Partition result;
  while ( unAssignedCluster.size() ) {  // as long as there is still object remains to be classified
    std::vector<CAST_Item> openCluster; // creates a new cluster
    resetAffinity( unAssignedCluster ); // reset the affinity related to the remaining objects
    bool changesOccurred = true; 
    while (changesOccurred && unAssignedCluster.size()) { // beging organizing objects
      changesOccurred = false;
      AffinityCompute maxCompute;
      AffinityCompute minCompute;
      while ( unAssignedCluster.size() ) { // tries to put the best-affinity objects in the open cluster
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
