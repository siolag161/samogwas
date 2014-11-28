 
#include "clustering/louvain/community.hpp"

#include "clustering/partition.hpp"
#include "clustering/louvain/graph.hpp"

#include <assert.h>     /* assert */

namespace samogwas
{
namespace louvain {

/**
 */
void Network::initialize() {
  m_labelSet.clear(); // clear labels
  m_index2Label.clear();

  community_member_counts.clear();
  community_member_counts.resize( nbrNodes(), 0); //

  in_weights.clear();
  in_weights.resize( nbrNodes(), 0);

  tot_linked_weights.clear();
  tot_linked_weights.resize( nbrNodes(), 0);
  for ( auto i = 0; i < nbrNodes(); ++i ) {
    addNode(i,i);
  }

  curr_modularity = INVALID_MODULARITY;
}

/**
 */
double Network::totalWeights() const {
  return graph->totalWeights();
}


/**
 */
CommunityIndex Network::communityOf( const NodeIndex& node ) const {
  // printf("network - communityOf: %d\n", node );
  return getLabel(node);
}

bool Network::sameCommunity( const NodeIndex& i, const NodeIndex& j ) const {
  return communityOf(i) == communityOf(j);
}


std::vector<NodeIndex> Network::membersOf( const CommunityIndex& comm ) const {
  assert( comm < nbrCommunities());
  std::vector<NodeIndex> members;

  for (int i = 0; i < nbrCommunities(); ++i) {
    if (communityOf(i) == comm ) {
      members.push_back(i);
    }
  }

  return members;
}

/**
 */
double Network::weight2community( const NodeIndex& node, const CommunityIndex& target ) const {
  
}

/**
 */
double Network::modularity() {

   if ( totalWeights() <= 0 ) {
     return INVALID_MODULARITY;
   }
  curr_modularity = 0.0;
  double tw2 = totalWeights()*2;

  for ( auto comm: communities() ) {
    if ( nbrCommunities() == 1) printf("sing fucking leton: %f\n", tot_linked_weights[comm]);
    curr_modularity += in_weights[comm] / tw2 - (tot_linked_weights[comm]/tw2)*(tot_linked_weights[comm]/tw2);
  }  
  return curr_modularity;
}


/**
 */
void Network::moveNode( const NodeIndex& node, const CommunityIndex& target,
                        const double old_shared_weights, const double new_shared_weights ) {
  CommunityIndex from = communityOf(node);
  if ( from == target ) return;
  removeNode(node, old_shared_weights);
  addNode(node, target, new_shared_weights);
}
 

/**
 */
void Network::moveNode( const NodeIndex& node, const CommunityIndex& target ) {
  CommunityIndex from = communityOf(node);
  if ( from == target ) return;
  double old_shared_weights = sharedWeights( node, from );
  double new_shared_weights = sharedWeights( node, target );
  moveNode( node, target, old_shared_weights, new_shared_weights);
}
  
double Network::modularityGain( const NodeIndex& node,
                                const CommunityIndex& comm,
                                const double shared_weights) const
{
  double gain = 0.0;
  if ( communityOf(node) == comm ) return gain;
  double sw = shared_weights;
  if ( shared_weights < 0 ) sw = sharedWeights(node,comm);
  double tw2 = 2*totalWeights();
  gain = sw - tot_linked_weights[comm]*linkedWeights(node)/tw2;

  // printf("gain of moving %d -> %d, %f. tot: %f, linked: %f\n", node, comm, gain,
  //        tot_linked_weights[comm], linkedWeights(node));
  return gain/totalWeights();
}

double Network::modularityLoss( const NodeIndex& node,
                                const double shared_weights ) const
{
  double loss = 0.0;
  const CommunityIndex comm = communityOf(node);
  double sw = shared_weights;
  if ( shared_weights < 0 ) sw = sharedWeights(node,comm);

    
  double tw2 = 2*totalWeights();
  if ( shared_weights < 0 ) sw = sharedWeights( node,comm );
  loss = -sw + (tot_linked_weights[comm] - linkedWeights(node) )*linkedWeights(node)/tw2;

  return loss/totalWeights();
}

/**
 */
size_t Network::nbrCommunities() const {
  return communities().size();
}

size_t Network::nbrNodes() const {
  return graph->nbrNodes();
}


size_t Network::nbrLinks() const {
  return graph->nbrLinks();
}


////////////////////////////////////
/** forces node --> community
 *
 */
void Network::setCommunity( const NodeIndex& node, const CommunityIndex& comm ) {
  setLabel(node, comm);
}



/**
 */
void Network::removeNode( const NodeIndex& node,
                          const double shared_weights) {
  CommunityIndex comm = communityOf(node);
  if ( comm >= 0 && comm < nbrNodes() ) {
    setCommunity(node, TEMP_COMMUNITY);
    community_member_counts[comm]--;
    if (community_member_counts[comm] == 0) {
      removeLabel(comm);
    }

    double tmp = tot_linked_weights[comm];
    tot_linked_weights[comm] -= linkedWeights(node);
    in_weights[comm] -=  2*shared_weights + graph->selfLoopWeight(node);
    // printf("removing %d from %d, old_tot: %f, new_tot: %f \n", node, comm, tmp, tot_linked_weights[comm]);
  }
}

/**
 */
void Network::addNode(const NodeIndex& node,
                      const CommunityIndex& comm,
                      const double shared_weights) {
  if ( comm >= 0 && comm < nbrNodes() ) {
    double tmp = in_weights[comm];
    setCommunity(node,comm);
    community_member_counts[comm]++;
    tot_linked_weights[comm] += linkedWeights(node);
    in_weights[comm] += 2*shared_weights + graph->selfLoopWeight(node);
  }  
}



double Network::sharedWeights( const NodeIndex& node, const CommunityIndex& comm ) const {
  double sw = 0.0;
  for ( auto vp = graph->adjacentNodes(node); vp.first != vp.second; ++vp.first ) {
    NodeIndex adj_node = *vp.first;
    if ( communityOf(adj_node) ==  comm ) {
      sw += graph->weight(node, adj_node);
    }
  }
  return sw;
}


double Network::interCommunityWeight( const CommunityIndex& i,
                                      const CommunityIndex& j ) {
  if ( i == j ) return in_weights[i];

  double inter_weight = 0.0;
  for ( auto vp = graph->allLinks(); vp.first != vp.second; ++vp.first ) {
    NodeIndex l = graph->source(*vp.first), r = graph->target(*vp.first);
    CommunityIndex c_l = communityOf(l), c_r = communityOf(r);
    if ( c_l == i && c_r == j || c_l == j && c_r == i ) { 
      inter_weight += graph->weight(l,r);
    }
  }
  return inter_weight;
}


} // samogwas
} 








