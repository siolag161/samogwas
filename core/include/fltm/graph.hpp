/****************************************************************************************
 * File: Graph.hpp
 * Description: This module provides the structure for the FLTM model.
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 13/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_GRAPH_HPP

#define SAMOGWAS_GRAPH_HPP

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include "pl.h"


namespace samogwas
{

typedef plJointDistribution JointDistribution;
typedef int IndexType;
typedef std::string LabelType;
typedef plSymbol RandomVariable;
typedef std::map<IndexType, LabelType> Index2Label;
typedef std::map<LabelType, IndexType> Label2Position;
typedef std::map<LabelType, IndexType> Label2Index;

/** Forward declaration for Node class as it is part of the Graph declaration
 *
 */
struct Node;

/* 
 *
 */
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS,
                                Node, boost::no_property> Graph;


/**
 * Node 
 */
struct Node {

  ~Node() { }
  
  IndexType index; 
  LabelType label; // cache for variable.name()
  
  bool isLeaf;
  int position; // physical position on the genome
  int level; // indicates the level to which this node belongs.

  Graph* graph; // reference to its graph

  RandomVariable variable; // represents the underlying random variable.
  JointDistribution jointDistribution; // joint distribution of the Naive Bayes Model rooted in the variable
  
  Node(): isLeaf(true), position(-1), level(0), graph(NULL) {}

  inline Node& setDistribution(JointDistribution &jointDistri);
  inline Node& setPosition(Label2Index &label2Index, plComputableObjectList &objectList);
  inline Node& setLevel(Label2Index &label2Index, plComputableObjectList &objectList);
  inline Node& setVariable(plComputableObjectList &objectList);
  
  /* Requires that the node is a latent variable.
   * Sets the following properties:
   *    - parent graph (from parameter)
   *    - joint distribution (from parameter)
   *    - variable (deduced from parameters)
   *    - position (deduced from parameters)
   *    - level (deduced from parameters).
   *
   * Requires Label2Index which maps the label of a variable to its index.
   */
  Node& setupProperties(Graph* graph,
                        JointDistribution& jointDist,
                        Label2Index& label2Index ) {
    plComputableObjectList objectList = jointDist.get_computable_object_list();

    setGraph(graph);
    setDistribution(jointDist);
    setVariable(objectList);
    setLabel(label2Index);
    setPosition(label2Index, objectList);
    setLevel(label2Index, objectList);
    return (*this);
    /**
    return setGraph(graph).setDistribution(jointDist)
            .setVariable(objectList).setLabel(label2Index)
            .setPosition(label2Index, objectList).setLevel(label2Index, objectList);
    */
  }
  

  Node& setGraph(Graph* graph) { this->graph = graph; return *this; }

  Node& setLabel(Label2Index& label2Index) {
    label = variable.name();
    label2Index[label] = this->index;
    return *this;
  }

  unsigned getIndex() const {
    return index;
  }
};

//////////////////////////////////////////////////////////////////
/**
 *
 */
typedef Graph::vertex_descriptor vertex_t;

/**
 *
 */
typedef Graph::edge_descriptor edge_t;

typedef boost::graph_traits< Graph >::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<Graph>::edge_iterator edge_iterator;


typedef std::pair<vertex_iterator, vertex_iterator> vertex_iter_pair;
typedef std::pair<edge_iterator, edge_iterator> edge_iter_pair;

/**
 *
 */
inline vertex_t createVertex( Graph& graph,
                              const unsigned &cardinality,
                              const bool isLeaf,                              
                              const std::string &label = "",
                              const unsigned &position = -1,
                              const unsigned &level = -1) {
  
  vertex_t vertexId = boost::add_vertex(graph); // adds a new Node to the graph and returns the newly added node's index.

  Node &newNode = graph[vertexId];
  newNode.variable = RandomVariable(label, plIntegerType(0, cardinality - 1));
  newNode.isLeaf = isLeaf;
  newNode.position = position; // physical position on the genome
  newNode.label = label;
  newNode.index = vertexId;
  newNode.graph = &graph;
  newNode.level = level;
  return vertexId;
}


//////////////////////////////////////////////////////////////////////////////////////////
Node& Node::setDistribution(JointDistribution &jointDist) {
  this->jointDistribution = jointDist;
  return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////
/**
 *
 */
Node& Node::setVariable(plComputableObjectList &objectList) {
  unsigned latentVarIdx = 0;
  this->variable = objectList[latentVarIdx].get_variables()[0];
  // objectList is the collection of distributions for
  // the Naive Bayes Model, augmented with annotations.
  // Example: P(X1,X2,Z) = P(Z)*P(X1|Z)*P(X2|Z).
  // ObjectList[0] = P(Z), ObjectList[1] = P(X1|Z), ObjectList[2] = P(X2|Z).
  // P(Z).get_variables() = {Z} ,P(X1|Z).get_variables() = {Z,X1}.
  
  return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////
// See above or a description of objectList.
// Computes the fictive position of the latent variable as the average of the children's positions.
Node& Node::setPosition(Label2Index &label2Index, plComputableObjectList &objectList) {
 
  unsigned nbrChildren = objectList.size() - 1;
  
  if (nbrChildren > 0) {    
    unsigned totalPos = 0;
    for (unsigned child = 0; child < nbrChildren; ++ child) {
      LabelType label = objectList[child].get_variables()[0].name();
      IndexType index = label2Index[label];
      totalPos += (*graph)[index].position; 
    }
    
    this->position = totalPos / nbrChildren;
  }

  return *this;
  
}

//////////////////////////////////////////////////////////////////////////////////////////
// See above or a description of objectList.
Node& Node::setLevel(Label2Index &label2Index, plComputableObjectList &objectList) {
  unsigned nbrChildren = objectList.size() - 1;

  int maxChildrenLevel = 0;
  for (unsigned child = 1; child <= nbrChildren; ++ child) {      
    LabelType label = objectList[child].get_variables()[0].name();
    IndexType index = label2Index[label];

    if (maxChildrenLevel < (*graph)[index].level) {
      maxChildrenLevel = (*graph)[index].level;
    }
  }
        
  this->level = maxChildrenLevel + 1;
  
  return *this;  
}


} // namespace samogwas ends here.

/****************************************************************************************/
#endif // SAMOGWAS_GRAPH_HPP
