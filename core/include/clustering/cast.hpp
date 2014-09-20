/****************************************************************************************
 * File: cast.hpp
 * Description: CAST (Clustering Affinity Search Technique) is a clustering algorithm proposed in Ben Dor et al. 1999
 * The idea of the algorithm is to represent the true underlying data as a probabilistic graph and a cluster in this graph is modeled as a clique in this structure. We introduce errors by randomly adding and removing edges in the graph betweeen pairs of nodes with a certain probability alpha. If the obtained cluster has a size of at least t*N for some constant t and N is the total number of items in the dataset then we may have found the cluster with a high confidence.
 * The implemention is based on the heuristic version described in the original paper where only one parameter is used. The idea is while we still have items that remain to be clustered (grouped together in a set called unassignedCluster), we try to generate a new cluster (openCluster) one by one. An item is added or removed from this openCluster if its relative affinity is greater or lower than the parameter t, respectively.
 *
 * @ref: Ben Dor et al. 1999, Ben-Dor, Amir, Ron Shamir, and Zohar Yakhini. "Clustering gene expression patterns." Journal of computational biology 6.3-4 (1999): 281-297.
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 11/07/2014                                                       

 ***************************************************************************************/
#ifndef SAMOGWAS_CAST_HPP // CS Homogeneize the names all throughout the project.
#define SAMOGWAS_CAST_HPP

#include "distance/similarity.hpp"

#include "clustering.hpp"
#include "partition.hpp"

namespace samogwas
{


/** The CAST methods returns a partition of the dataset, much like other type of AlgoClust. Under the hood, during computation make use of the private
 * class CAST_Item to represent an item in the cluster. CAST_Item is not supposed to be used by outside class.
 * CS CS CLARIFY in the form of non-overlapping set sets ??? of <objects of type> CAST_Item ???.
 *  A CAST_Item is a structure which holds an index of the object ( <CS suggestion index relative to the> relative to the original dataset ) and
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
template<typename SimiMatrix> // CS I am lost 
                              // I do not know if your SimiMatrix is elsewhere your compMatrix.
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
    for ( int i = 0; i < this->compMatrix->size(); ++i ) { // CS Nobody can guess that comp is a matrix. compMat should be used.
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
// CS What about the separation rules? Where to insert /////////// or /*********/?

struct AffinityCompute { // CS heterogeneous notation: resetAffinity versus AffinityCompute @pdt: not the same, one is a class, the other is a method
  template<typename Compare>
  const int operator()(const std::vector<CAST_Item>& clust, Compare comp) const;    
};

// CS What about the separation rules? Where to insert /////////// or /*********/?

/** Resets the current affinity of the cluster given as parameter
 *
 */
inline void resetAffinity( std::vector<CAST_Item>& cluster);
// What is the link with the above inline void resetAffinity?
//
// CS I do not see that there is such a step in the original CAST algorithm. @pdt: indeed, everytime a new openCluster is considered, gotta reset all the affinity relative to this one.
// Instead I see : When a new cluster Copen is started, the initial affinities of all genes are 0 since Copen is empty.
// CS Is it this step you describe here?
// http://faculty.washington.edu/kayee/cluster/algs.pdf
// Is it resetAffinity for the clusterofUnassignedObjects?


// CS is unassignedCluster the equivalent of COpen? // pdt: not the same. COpen is the new open cluster to be formed, unassigned is a group of remaining objects
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


/** Implements the heuristic CAST algorithm following the original paper mentioned in the filehead.
 *  It requires a similarity matrix ( or any method that returns a similarity between a pair of item index)
 *  as input. It also requires a parameter t as the affinity threshold which determines whether we removes
 *  or add an item to the current candidate cluster. This parameter t also determines when a cluster is stabilizes.
 */ 
template<typename SimiMatrix>
Partition CAST<SimiMatrix>::run( std::vector<CAST_Item>& unAssignedCluster ) {
  Partition result; // CS growingPartition ?
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
                                           // CS successively tries to assign the objects with the best affinities to openCluster
        int maxAffIdx = maxCompute( unAssignedCluster, std::greater<double>() );
        if ( unAssignedCluster.at(maxAffIdx).affinity >= thresCAST * openCluster.size() ) {
          changesOccurred = true;
          addGoodItem( unAssignedCluster, openCluster, *this->compMatrix, maxAffIdx );
        } else {
          break;
        }
      }
      while( unAssignedCluster.size() ) { // tries to remove the worst-affinity objects from the open cluster
                                          // CS successively tries to remove the objects with the worst affinities from openCluster
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
      result.cluster( item.index, cluster_id ); // CS The procedure's name is ambiguous.
    }
  }  
  return result;  
}

// sets all the affinity values to zeroes
template<typename SimiMatrix>
void CAST<SimiMatrix>::resetAffinity( std::vector<CAST_Item>& cluster ) { // CS resetAffinityValuesToZero
  for ( std::vector<CAST_Item>::iterator it = cluster.begin(); it != cluster.end(); ++it ) {
    it->affinity = 0.0;
  } 
}

////////////////////////////////////////////////////////////////////////////////////////
// CS What is the rule for using //////////// or /*****************/?
/** Searches for the object in given cluster which has the optimal (best/worst) affinity
 *
 */
template<typename Compare>
const int AffinityCompute::operator()( const std::vector<CAST_Item>& clust, Compare comp) const
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

////////////////////////////////////////////////////////////////////////////////////////

/// Adds a potentially good candidate item to openCluster
void addGoodItem( std::vector<CAST_Item>& unAssignedCluster, std::vector<CAST_Item>& openCluster, 
                  CompMatrix& simCompute, const int clusterIdx ) {  
  updateClustersAffinity( unAssignedCluster, openCluster, 
                          simCompute, clusterIdx );  
}

/// removes a bad item from the cluster CS from openCluster ?
void removeBadItem( std::vector<CAST_Item>& unAssignedCluster, std::vector<CAST_Item>& openCluster, 
                    CompMatrix& simCompute, const int clusterIdx ) {
  updateClustersAffinity( openCluster, unAssignedCluster, 
                          simCompute, clusterIdx); 
}

// updates the affinity after changes are made to reflect current assignment
// CS updates the affinity values in sourceCluster and targetCluster after the move of indeXItem
void updateClustersAffinity( std::vector<CAST_Item>& sourceCluster, std::vector<CAST_Item>& targetCluster,
                             CompMatrix& simCompute, 
                             const int clusterIndex ) // CS bad denomination
                                                      // indexOfItemInSourceCluster
{
  const CAST_Item item = sourceCluster.at(clusterIndex);
  // CS not optimized?
  // CS const int itemIndex = sourceCluster.at(clusterIndex).index;
  moveItemBetweenClusters( sourceCluster, targetCluster, clusterIndex ); 
  for (int i = 0; i < sourceCluster.size(); i++) {
    sourceCluster[i].affinity += simCompute( sourceCluster[i].index, item.index ); // CS I am lost
                                                                                   // Is type CompMatrix for ComparisonMatrix
                                                                                   // or is it for simComputationFromMatrix?
                                                                                   // suggestion for type identifier: CompSimFromMatrix?
                                                                                   // + itemIndex instead of item.index
  }

  for (int i = 0; i < targetCluster.size(); i++) {
    targetCluster[i].affinity += simCompute( targetCluster[i].index, item.index ); // CS itemIndex 
  } 

}

/// Moves item from source cluster to target cluster
void moveItemBetweenClusters( std::vector<CAST_Item>& source,
                              std::vector<CAST_Item>& target,
                              const int clusterIndex ) { // indexOfItemInSourceCluster
  target.push_back( source.at(clusterIndex) );
  //source.remove(clusterIndex);
  source.erase( source.begin() + clusterIndex ); // CS I do not understand this argument
                                                 // CS C++ vector erase Removes the element at pos.
                                                  // CS I would have expected indexOfItemInSourceCluster
}


}

#endif // SAMOGWAS_CAST_HPP
