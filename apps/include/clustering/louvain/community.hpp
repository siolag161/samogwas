/****************************************************************************************
 * File: community.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 21/11/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_LOUVAIN_COMMUNITY_HPP
#define SAMOGWAS_LOUVAIN_COMMUNITY_HPP

#include <set>
#include "clustering/partition.hpp"
#include "graph.hpp"
namespace samogwas
{

namespace louvain {

//typedef Index NodeIndex;
typedef Partition::Label CommunityIndex;

class Network: public Partition {
  
 public:
  static const int TEMP_COMMUNITY = -1;
  static const int INVALID_MODULARITY = -2;
  
  /**
   */
  Network(std::shared_ptr<Graph> g): graph(g) { initialize(); }

  void setGraph(std::shared_ptr<Graph> g) { graph = g; initialize(); }

  /**
   */
  void removeNode(const NodeIndex& node, const double shared_weights = 0);

  void moveNode( const NodeIndex& node, const CommunityIndex& target,
                          const double old_shared_weights, const double new_shared_weights );
  
      /**
       */
  void moveNode(const NodeIndex& node, const CommunityIndex& target);

  
  /**
   */
  void addNode(const NodeIndex& node, const CommunityIndex& comm, const double shared_weights = 0);

  /**
   */
  double modularityLoss( const NodeIndex& node,
                         // const double own_shared, 
                         const double shared_weights = -1.0) const;

  /**
   */
  double modularityGain( const NodeIndex& node,
                         const CommunityIndex& target,
                         // const double own_shared, 
                         const double shared_weights = -1.0) const;

  // /**
  //  */
  // double modularityGain( const NodeIndex& node,
  //                        const CommunityIndex& comm,
  //                        const double shared_weights = -1.0) const;
  /**
   */
  CommunityIndex communityOf( const NodeIndex& node ) const;


  
  Graph::AdjIteRng adjacentNodes( const NodeIndex& i ) const { return graph->adjacentNodes(i); }

  /**
   */
  bool sameCommunity( const NodeIndex& i, const NodeIndex& j ) const;
      
  /**
   */
  double modularity();

  /**
   */
  size_t nbrCommunities() const;

  /**
   */
  size_t nbrNodes() const;

  
  /**
   */
  size_t nbrLinks() const;
  
  //protected:

  std::vector<NodeIndex> membersOf( const CommunityIndex& comm ) const;
  
  /**
   */
  void initialize();


  /** sum of weights from a node to all the nodes a given community
   */
  double weight2community( const NodeIndex& node, const CommunityIndex& comm ) const;
  
  /**
   */
  double totalWeights() const;


  // /**  given a node, computes the total weights of other nodes connected to this one
  //  *
  //  */
  double linkedWeights( const NodeIndex& node ) const { return graph->linkedWeights(node); }

  //
  void setCommunity( const NodeIndex& node, const CommunityIndex& comm );

  double interCommunityWeight( const CommunityIndex& i, const CommunityIndex& j);

  
  virtual size_t nbrClusters() const { return m_labelSet.size(); }
  virtual size_t nbrItems() const { return graph->nbrNodes(); }
  virtual void setLabel( NodeIndex node, CommunityIndex comm ) {  //@doto: setLabel
    if (comm != TEMP_COMMUNITY) m_labelSet.insert(comm);
    m_index2Label[node] = comm;
  }

  const std::set<CommunityIndex>& communities() const {
    return m_labelSet;
  }
  
  double sharedWeights( const NodeIndex& node, const CommunityIndex& comm ) const;
  

 public:
  std::shared_ptr<Graph> graph;

  // int nbr_communities;
  std::vector<int> community_member_counts; // keeps track of the nbr of member per comunity
  
  // modularities & links
  std::vector<double> in_weights; // total weight of links (including self-loops) inside each community
  std::vector<double> tot_linked_weights;


  double curr_modularity;

  // double valid_values;
};

} // louvain ends here.

} // namespace samogwasends here. samogwas

/****************************** IMLEMENTATION BELOW THIS POINT **************************/
namespace samogwas
{


} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_LOUVAIN_COMMUNITY_HPP
