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
  const int TEMP_COMMUNITY = -1;
  const int INVALID_MODULARITY = -2;
  typedef double Weight;
  /**
   */
  Network(std::shared_ptr<Graph> g): graph(g) { initialize(); }

  void setGraph(std::shared_ptr<Graph> g) { graph = g; initialize(); }

  /**
   */
  void removeNode(const NodeIndex& node, const Weight shared_weights = 0);

  void moveNode( const NodeIndex& node, const CommunityIndex& target,
                          const Weight old_shared_weights, const Weight new_shared_weights );
  
      /**
       */
  void moveNode(const NodeIndex& node, const CommunityIndex& target);

  
  /**
   */
  void addNode(const NodeIndex& node, const CommunityIndex& comm, const Weight shared_weights = 0);

  /**
   */
  Weight modularityLoss( const NodeIndex& node,
                         // const Weight own_shared, 
                         const Weight shared_weights = -1.0) const;

  /**
   */
  Weight modularityGain( const NodeIndex& node,
                         const CommunityIndex& target,
                         // const Weight own_shared, 
                         const Weight shared_weights = -1.0) const;

  // /**
  //  */
  // Weight modularityGain( const NodeIndex& node,
  //                        const CommunityIndex& comm,
  //                        const Weight shared_weights = -1.0) const;
  /**
   */
  CommunityIndex getCommunity( const NodeIndex& node ) const;

  Graph::AdjIteRng adjacentNodes( const NodeIndex& i ) const { return graph->adjacentNodes(i); }


  Graph::LinkIteRng allLinks() const { return graph->allLinks(); }
    
  /**
   */
  NodeIndex source( const Graph::LinkIndex& link ) const { return graph->source(link); }

  /**
   */
  NodeIndex target( const Graph::LinkIndex& link ) const { return graph->target(link); }
  
  /**
   */
  bool sameCommunity( const NodeIndex& i, const NodeIndex& j ) const;
      
  /**
   */
  Weight modularity();

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


  // /** sum of weights from a node to all the nodes a given community
  //  */
  // Weight weight2community( const NodeIndex& node, const CommunityIndex& comm ) const;
  
  /**
   */
  Weight totalWeights() const;

    /**
   *
   */
  Weight weight( const NodeIndex& i, const NodeIndex& j) const {
    return graph->weight(i,j);
  }

  /**
   *
   */
  Weight weight( const Graph::LinkIndex& j) const { return graph->weight(j); }


  // /**  given a node, computes the total weights of other nodes connected to this one
  //  *
  //  */
  Weight linkedWeights( const NodeIndex& node ) const { return graph->linkedWeights(node); }

  //
  void setCommunity( const NodeIndex& node, const CommunityIndex& comm );

  Weight interCommunityWeight( const CommunityIndex& i, const CommunityIndex& j);

  Weight innerCommunityWeight( const CommunityIndex& i );
  
  virtual size_t nbrClusters() const { return m_labelSet.size(); }
  virtual size_t nbrItems() const { return graph->nbrNodes(); }
  virtual void setLabel( NodeIndex node, CommunityIndex comm ) {  //@doto: setLabel
    if (comm != TEMP_COMMUNITY) m_labelSet.insert(comm);
    m_index2Label[node] = comm;
  }

  const std::set<CommunityIndex>& communities() const {
    return m_labelSet;
  }
  
  Weight sharedWeights( const NodeIndex& node, const CommunityIndex& comm ) const;
  

  // void buildDendogram(const Graph& curr_graph);
 public:
  std::shared_ptr<Graph> graph;

  // int nbr_communities;
  std::vector<int> community_member_counts; // keeps track of the nbr of member per comunity
  
  // modularities & links
  std::vector<Weight> in_weights; // total weight of links (including self-loops) inside each community
  std::vector<Weight> tot_linked_weights;


  Weight curr_modularity;

  // Weight valid_values;
};

} // louvain ends here.

} // namespace samogwasends here. samogwas

/****************************** IMLEMENTATION BELOW THIS POINT **************************/
namespace samogwas
{


} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_LOUVAIN_COMMUNITY_HPP
