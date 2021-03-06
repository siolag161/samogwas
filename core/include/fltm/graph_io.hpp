/****************************************************************************************
 * File: GraphIO.hpp lo
 * Description: This module provides methods to save an FLTM Graph into one or two formatted files during the
 * ------------ FLTM construction process; an additional file relative to positions is output in every case.
 * -----------  This module also provides methods to load an FLTM graph from such
 * -----------  formatted files output by the FLTM construction process.
 * -----------  Three formats are available (SINGLE, TULIP, BN).
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 03/04/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_GRAPH_GRAPHIO_HPP
#define SAMOGWAS_GRAPH_GRAPHIO_HPP

#include "graph.hpp"
#include "utils/csv_parser.hpp"
#include <map>
#include <string>
#include <pl.h>
#include <iomanip>
#include <boost/lexical_cast.hpp>

namespace samogwas
{

static const char SEPARATOR = ',';

enum FULL { FULL_ID = 0, FULL_LABEL, FULL_LATENT, FULL_PARENT, FULL_LEVEL, FULL_POSITION, FULL_CARDINALITY };
enum LAB_POS { LP_ID = 0, LP_LABEL, LP_POSITION };

// format SINGLE (one file, structure)
enum SINGLE { SINGLE_ID = 0, SINGLE_LATENT, SINGLE_PARENT, SINGLE_LEVEL, SINGLE_CARDINALITY };
// by convention, SINGLE_PARENT set to -1 denotes a root.

// format TULIP (two files, structure)
enum TULIP_VERTICES { TULIP_ID = 0, TULIP_LATENT, TULIP_LEVEL, TULIP_CARDINALITY };
enum TULIP_EDGES { /*TULIP_ID = 0,*/ TULIP_PARENT_ID = 1 };
// by convention, TULIP_PARENT_ID set to -1 denotes a root.


// format BN (two files, structure and parameters)
enum BN_VERTICES { BN_LATENT_ID = 0, NBR_CHILDREN };

static const std::string ID = "id"; static const std::string LABEL = "label";
static const std::string LATENT = "latent"; static const std::string PARENT = "parent";
static const std::string LEVEL = "level"; static const std::string POSITION = "position";
static const std::string CARDINALITY = "cardinality";
static const std::string PARENT_ID = "parent_id";

typedef std::map<int, std::pair<std::string, int> > LabPosMap;

/*
 *  
 */
struct FLTMGraphReader {
  /** Loads the nodes' idendities, labels and positions from a formatted file and returns a
   * map that maps the identity to a pair of label and position.
   */
  inline LabPosMap readLabPos( const std::string labPosFileName ) const;  
};

/** 
 *
 */
struct SingleGraphLoad: public FLTMGraphReader {

    /** Reconstructs an FLTM graph from the SINGLE formatted file.
    */
  inline void operator()( Graph& graph, const std::string labPosFileName, const std::string singleFileName ) const;

  Graph operator()( const std::string labPosFileName, const std::string singleFileName ) const {
    Graph graph; (*this)(graph, labPosFileName, singleFileName); return graph;
  }
};

/** Loads a graph which represents the topology of the FLTM Model from two separate files: one containing information regarding edges
 *  and the other concerning the vertices.
 */
struct TulipGraphLoad: public FLTMGraphReader {

  inline void operator()( Graph& graph, const std::string labPosFileName,
                          const std::string vertexFileName, const std::string edgeFileName ) const;

  Graph operator()(  const std::string labPosFileName, const std::string vertexFileName, const std::string edgeFileName) const {
    Graph graph; (*this)( graph, labPosFileName, vertexFileName, edgeFileName ); return graph;
  }
};


/** Loads the Bayian network corresponding to the FLTM model from two files:
 *  - vertexFileName: information regarding vertices of the following (SINGLE) format: ID, LATENT, PARENT, LEVEL, CARDINALITY
 *  - distributionFileName:
 *     == example: Z has 2 children X1 and X2; card(Z) = 3; card(X1) = card(X2) = 2;
 *         Z_id 2 // 2 denotes the number of children.
 *         p(Z=0) p(Z=1) p(Z=2)
 *         X1_id
 *         P(x1=0 | Z=0) p(X1=1 | Z=0)
 *         p(x1 | Z=1)
 *         p(x1 | Z=2)
 *         X2_id
 *         P(x2=0 | Z=0) p(X2=1 | Z=0)
 *         p(x2 | Z=1)
 *         p(x2 | Z=2)
 *
 *         Y_id 7
 *         etc.
 */
struct BayesGraphLoad: public FLTMGraphReader {  
  inline void operator()( Graph& graph, const std::string labPosFileName,
                          const std::string vertexFileName, const std::string distributionFileName ) const;
  Graph operator()( const std::string labPosFileName, const std::string vertexFileName, const std::string distributionFileName ) const  {
    Graph graph;  (*this)( graph, labPosFileName, vertexFileName, distributionFileName );  return graph;
  }
};
///////////////////////////////////////////////////////////////

/** @todo: check usefulness
 *
 */
struct FLTMGraphWriter {  
  inline void writeLabPos( const Graph& graph, const std::string labPosFileName );  
};


/** Saves the graph which represents the topology of the FLTM Model into two files:
 *    - one containing the label-position information
 *    - the other containing: id, latent, id of the parent, level and cardinality. @todo: check FLTMGraphReader
 */
struct SingleGraphSave {
  inline void operator()( Graph& graph, const std::string singleFileName ) const;
};

/** Saves an FLTM Graph into two files:
  *    - one containing information regarding vertices with the following format: ID, LATENT, LEVEL
  *    - the other concerning the edges: ID, ID_PARENT.
  */
struct TulipGraphSave {  
  inline void operator()( const Graph& graph, const std::string vertexFileName, const std::string edgeFileName) const;
};


/** Saves the Bayesian network corresponding to the FLTM model into two files of the following formats:
 *  - vertexOutputFilename: ID, LABEL, LATENT, PARENT, LEVEL, POSITION, CARDINALITY
 *  - distributionFileName: see above.
 */
struct BayesGraphSave {  
  inline void operator()( const Graph& graph, const std::string vertexFileName, const std::string distributionFileName ) const;
};

} // namespace samogwas_graph ends here.

/****************************** IMPLEMENTATION BELOW THIS POINT **************************/
namespace samogwas
{

 // Format LP_ID = 0, LP_LABEL, LP_POSITION
LabPosMap FLTMGraphReader::readLabPos( const std::string labPosFileName ) const {
  LabPosMap lpMap;
  std::ifstream labPosFile(labPosFileName.c_str());    
  utility::CSVIterator<std::string> labPosLine(labPosFile);

  size_t id = 0;
  for( ; labPosLine != utility::CSVIterator<std::string>(); ++labPosLine ) {
    int position; std::string label;
    if (labPosLine->size() == 4) {
      // id = boost::lexical_cast<size_t>( (*labPosLine)[LP_ID] );
      label = (*labPosLine)[LP_LABEL] ;
      position = boost::lexical_cast<int>( (*labPosLine)[LP_POSITION] );
    } else {
      // id = boost::lexical_cast<size_t>( (*labPosLine)[LP_ID+1] );
      label = (*labPosLine)[LP_LABEL+1] ;
      position = boost::lexical_cast<int>( (*labPosLine)[LP_POSITION+1] );
    }
    lpMap[id] = std::pair<std::string, int>(label, position);
    id++;
  }
  return lpMap;
}

//////////////////////////////////////////////////////////////////////////////////////////
// SINGLE format: SINGLE_ID = 0, SINGLE_LATENT, SINGLE_PARENT, SINGLE_LEVEL, SINGLE_CARDINALITY
void SingleGraphLoad::operator()( Graph& graph, const std::string labPosFileName, const std::string singleFileName ) const {
  
  std::ifstream singleFile(singleFileName.c_str());  
  if (!singleFile) {
    return ;
  }

  LabPosMap labPosMap = readLabPos(labPosFileName);
    
  typedef std::map<int,int> EdgeMap;
  EdgeMap parentOf;
  utility::CSVIterator<std::string> singleLine(singleFile); ++singleLine;  // skips the header.
  for( ; singleLine != utility::CSVIterator<std::string>(); ++singleLine ) {
    size_t id = boost::lexical_cast<size_t>( (*singleLine)[SINGLE_ID] );   
    bool isLeaf = !(boost::lexical_cast<bool>( (*singleLine)[SINGLE_LATENT])) ;
    size_t parentId = boost::lexical_cast<size_t>( (*singleLine)[SINGLE_PARENT] );
    int level =  boost::lexical_cast<size_t>( (*singleLine)[SINGLE_LEVEL] );  
    size_t cardinality = boost::lexical_cast<size_t>( (*singleLine)[SINGLE_CARDINALITY] );      
    parentOf[id] = parentId;
    
    std::string label = labPosMap[id].first;
    int position = labPosMap[id].second;

    createVertex( graph, cardinality, isLeaf, label, position, level );
  }
  
  for ( EdgeMap::const_iterator it = parentOf.begin(); it != parentOf.end(); ++it ) {
    if (it->second != -1) { // if the node is not a root (By convention, -1 denotes that the node is a root.)
      boost::add_edge( it->second, it->first, graph );
    }
  }

  singleFile.close();
}

///////////////////////////////////////////////////////////////////////////////////////
// Tulip Vertex file's format: TULIP_ID, TULIP_LATENT, TULIP_LEVEL, TULIP_CARDINALITY
// Tulip edge file's format: TULIP_ID, TULIP_PARENT_ID
void TulipGraphLoad::operator()( Graph& graph,
                                 const std::string labPosFileName,
                                 const std::string vertexFileName,
                                 const std::string edgeFileName ) const {
  
  LabPosMap labPosMap = readLabPos(labPosFileName);
  
  std::ifstream vertexFile(vertexFileName.c_str());  
  std::ifstream edgeFile(edgeFileName.c_str());
  if (!edgeFile || !vertexFile) {
    return ;
  }

  utility::CSVIterator<std::string> vertexLine(vertexFile); ++vertexLine;  // skips the header.
  for( ; vertexLine != utility::CSVIterator<std::string>(); ++vertexLine ) {         
    size_t id = boost::lexical_cast<size_t>( (*vertexLine)[TULIP_ID] );   
    bool isLeaf = !(boost::lexical_cast<bool>( (*vertexLine)[TULIP_LATENT])) ;
    int level =  boost::lexical_cast<size_t>( (*vertexLine)[TULIP_LEVEL] );  
    size_t cardinality = boost::lexical_cast<size_t>( (*vertexLine)[TULIP_CARDINALITY] );      
    
    std::string label = labPosMap[id].first;
    int position = labPosMap[id].second;

    createVertex( graph, cardinality, isLeaf, label, position, level );
  }
  
  
 utility::CSVIterator<std::string> edgeLine(edgeFile); ++edgeLine; // skips the header.
  for(; edgeLine != utility::CSVIterator<std::string>(); ++edgeLine ) {    
     
    int sourceId = boost::lexical_cast<int>( (*edgeLine)[TULIP_PARENT_ID] );
    int targetId = boost::lexical_cast<int>( (*edgeLine)[TULIP_ID] );
    boost::add_edge( sourceId, targetId, graph );
  }
  edgeFile.close();
}

///////////////////////////////////////////////////////////////////////////////////////

/** BN vertex file: see SINGLE format.
 *  BN distribution file: see above.
 */
void BayesGraphLoad::operator()( Graph& graph,
                                 const std::string labPosFileName,
                                 const std::string vertexFileName,
                                 const std::string distributionFileName ) const {  

  std::cout << "begin loading label..." << std::endl;
  LabPosMap labPosMap = readLabPos(labPosFileName);
  std::cout << "end loading label..." << std::endl;

  std::ifstream vertexFile(vertexFileName.c_str()), distributionFile(distributionFileName.c_str());

  utility::CSVIterator<std::string> vertexLine(vertexFile); ++vertexLine; // skips header.
  std::cout << "begin loading vertices\n";
  for( ; vertexLine != utility::CSVIterator<std::string>(); ++vertexLine ) {
    size_t id = boost::lexical_cast<size_t>( (*vertexLine)[TULIP_ID] );
    bool isLeaf = !(boost::lexical_cast<bool>( (*vertexLine)[TULIP_LATENT]));
    int level =  boost::lexical_cast<size_t>( (*vertexLine)[TULIP_LEVEL] );
    size_t cardinality = boost::lexical_cast<size_t>( (*vertexLine)[TULIP_CARDINALITY] );
    std::string label = labPosMap[id].first;
    int position = labPosMap[id].second;
    createVertex( graph, cardinality, isLeaf, label, position, level );
  }
  std::cout << "end loading vertices: " << boost::num_vertices(graph) <<  "\n\n";
  std::cout << "begin loading dist\n\n";

  utility::CSVIterator<std::string> distributionLine(distributionFile);
  while ( distributionLine != utility::CSVIterator<std::string>() ) {
    plVariablesConjunction variables; // holds child variables
    plComputableObjectList jointDistri;
    size_t latentId =  boost::lexical_cast<size_t>( (*distributionLine)[BN_LATENT_ID] );
    size_t nbrChildren = boost::lexical_cast<size_t>( (*distributionLine)[NBR_CHILDREN] );
    Node& latentNode = graph[ latentId ]; ++distributionLine; // reads next line.

    std::vector< plProbValue > probValuesZ;
    for ( size_t latentVal = 0; latentVal < latentNode.variable.cardinality(); ++latentVal) { // loads probability table
                                                                                              // for the latent var
      probValuesZ.push_back( boost::lexical_cast<plProbValue>( (*distributionLine)[latentVal] ) );
    }
    const plProbTable probTabZ(latentNode.variable, probValuesZ); ++distributionLine;    
    for ( size_t child = 0; child < nbrChildren; ++child ) {
      size_t childId = boost::lexical_cast<size_t>( (*distributionLine)[BN_LATENT_ID] ); ++distributionLine;
      Node& childNode = graph[ childId ]; variables ^= childNode.variable;     
      plDistributionTable distTab_Xi_Z ( childNode.variable, latentNode.variable );
      for ( size_t latentVal = 0; latentVal < latentNode.variable.cardinality(); ++latentVal ) {
        std::vector< plProbValue > probValuesXiZ_vals;
        for ( size_t childVal = 0; childVal < childNode.variable.cardinality(); ++childVal ) {
          probValuesXiZ_vals.push_back( boost::lexical_cast<plProbValue>( (*distributionLine)[childVal] ) );
        }
        distTab_Xi_Z.push( plProbTable( childNode.variable, probValuesXiZ_vals), (int)latentVal );
        ++distributionLine;
      }
      jointDistri *= distTab_Xi_Z; // adds the conditional table to result
      boost::add_edge( latentId, childId, graph );
    }
    ++distributionLine;
    latentNode.jointDistribution = plJointDistribution(latentNode.variable ^ variables, probTabZ * jointDistri);
  }
  distributionFile.close(); vertexFile.close();  
}

////////////////////////////////// WRITER BELOW THIS //////////////////////////////////////

/*void FLTMGraphWriter::writeLabPos( const Graph& graph,
                                   const std::string labPosFileName ) {
  std::ofstream labPosFile(labPosFileName.c_str());  
  vertex_iterator vi, vi_end;
  labPosFile << ID << SEPARATOR << LABEL << SEPARATOR << POSITION << SEPARATOR << std::endl;  // writes header
  for ( boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi ) {
    int vertex = *vi;
    labPosFile << graph[vertex].index << SEPARATOR
               << graph[vertex].label << SEPARATOR
               << graph[vertex].position << std::endl; 
  }  
  labPosFile.close();
}*/

////////////////////////////////////////////////////////////////////////////////////////////
// SINGLE format: SINGLE_ID = 0, SINGLE_LATENT, SINGLE_PARENT, SINGLE_LEVEL, SINGLE_CARDINALITY
void SingleGraphSave::operator()( Graph& graph, const std::string singleFileName ) const {
  
  std::ofstream singleFile(singleFileName.c_str());
  
  // std::map<int, int> parent;
  // edge_iterator ei, ei_end;
  // for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei) {
  //   edgeMap[boost::target(*ei, graph)] = boost::source(*ei, graph);
  // }
  std::cout << "wrting graph of " << boost::num_edges(graph) << " edges " << std::endl;
  std::vector<int> parent( boost::num_vertices(graph), -1);
  for ( auto ep = boost::edges(graph); ep.first != ep.second; ++ep.first ) {
    int s = boost::source(*ep.first, graph), t = boost::target(*ep.first, graph);
    parent[t] = s;
  }
   
  vertex_iterator vi, vi_end;
  singleFile << ID << SEPARATOR <<  LATENT << SEPARATOR << PARENT
             << SEPARATOR << LEVEL << SEPARATOR << CARDINALITY << "\n";  // writes header
  for ( boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi ) {
    int vertex = *vi;
    singleFile << graph[vertex].index << SEPARATOR
               << !graph[vertex].isLeaf << SEPARATOR    
               << parent[vertex] << SEPARATOR
               << graph[vertex].level << SEPARATOR
               << graph[vertex].variable.cardinality() << std::endl;
    // printf("we have this %d->%d\n", parent[vertex], vertex);

  }  
  
  singleFile.close();
}


/////////////////////////////////////////////////////////////////////////////////////////
// Tulip Vertex file's format: TULIP_ID, TULIP_LATENT, TULIP_LEVEL, TULIP_CARDINALITY
// Tulip edge file's format: TULIP_ID, TULIP_PARENT_ID
void TulipGraphSave::operator()( const Graph& graph,
                                 const std::string vertexFileName,
                                 const std::string edgeFileName ) const {
 
  std::ofstream vertexFile(vertexFileName.c_str());
  vertex_iterator vi, vi_end;
  vertexFile << ID << SEPARATOR << LATENT << SEPARATOR << LEVEL << SEPARATOR << CARDINALITY << "\n";  // writes header
  for ( boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi ) {
    int vertex = *vi;
    vertexFile << graph[vertex].index << SEPARATOR
               << !(graph[vertex].isLeaf) << SEPARATOR
               << graph[vertex].level << SEPARATOR
               << graph[vertex].variable.cardinality()
               << std::endl;
  }  
  vertexFile.close();

  std::ofstream edgeFile(edgeFileName.c_str());
  edgeFile << ID << SEPARATOR << PARENT_ID << "\n"; // writes header
  edge_iterator ei, ei_end;
  for( boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei ) {
    edgeFile << boost::target(*ei, graph) << SEPARATOR << boost::source(*ei, graph) << std::endl;
  }
  edgeFile.close(); 
 
}

////////////////////////////////////////////////////////////////////////////////
/** vertex file: see SINGLE format.
*  distribution file: see above.
*/
void BayesGraphSave::operator()( const Graph& graph,
                                 const std::string vertexFileName,
                                 const std::string distributionFileName ) const {
  
  std::ofstream distributionFile(distributionFileName.c_str()), vertexFile(vertexFileName.c_str());
  vertex_iterator vi, vi_end;
  Label2Index label2Index;
  vertexFile << ID << SEPARATOR << LATENT << SEPARATOR << LEVEL << SEPARATOR << CARDINALITY << "\n";  // writes header
  std::cout << "saving vertices...\n";
  for ( boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi ) {
    int vertex = *vi;
    vertexFile << graph[vertex].index << SEPARATOR << !(graph[vertex].isLeaf) << SEPARATOR
               << graph[vertex].level << SEPARATOR << graph[vertex].variable.cardinality() << std::endl;
    label2Index[ graph[*vi].label ] = graph[*vi].index;
  }
  vertexFile.close();
  std::cout << "saving joint distribution...\n";
   
  for ( boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi ) {
    const Node& node = graph[*vi];                 
    if ( !node.isLeaf ) {
      plJointDistribution distribution = node.jointDistribution;
      plVariablesConjunction all_variables = distribution.get_variables(); // all variables (latent variable and its children)

      plVariablesConjunction variables; // child variables
      for (size_t i = 1; i <  all_variables.size(); ++i)
        variables ^=  all_variables[i]; // initializes child conjunction.

      plSymbol latentVar =  all_variables[0]; // latent variable
      distributionFile << node.index << SEPARATOR <<  variables.size() << std::endl;

      plComputableObjectList objLists = distribution.get_computable_object_list();

      plComputableObject probTableZ = objLists.get_distribution_over(latentVar); // distribution table for the latent variable
      int val;
      for ( val = 0; val < latentVar.cardinality() - 1 ; ++val ) {      
        distributionFile << std::fixed << std::setprecision(30) << probTableZ.compute( plValues().add(latentVar, val) ) << SEPARATOR; // P(latentVar = val)
      }

      distributionFile << std::fixed << std::setprecision(30) << probTableZ.compute( plValues().add(latentVar, val) )  << std::endl; // writes last probability value
      for ( size_t i = 0; i < variables.size(); ++i ) {
        plSymbol varX = variables[ i ]; // retrieves the child variable
        distributionFile << label2Index[ varX.name() ]  << std::endl; // writes child variable's id.
        plComputableObject compTableXZ = objLists.get_distribution_over( varX ); // conditional distribution P(X_i | Z)
        plDistributionTable& distTableXZ = static_cast<plDistributionTable&>( compTableXZ ); // casting P(X_i | Z) to derived class
        for ( val = 0; val < latentVar.cardinality(); ++val ) {
          int childVal;
          for ( childVal = 0; childVal < varX.cardinality() - 1; ++childVal ) { // for each value x of the child variable
            distributionFile << std::fixed << std::setprecision(30) << distTableXZ.compute( plValues().add(latentVar, val).add(varX, childVal) ) << SEPARATOR; // p(X_i = childVal | Z = val)
          }
          distributionFile << std::fixed << std::setprecision(30) << distTableXZ.compute( plValues().add(latentVar, val).add(varX, childVal) ) << std::endl;
        }
      }
      distributionFile << std::endl; // breaks the line, moves to the next latent variable.
    }
  }
 
  distributionFile.close();
}


} // namespace samogwas ends here.
/****************************************************************************************/
#endif // SAMOGWAS_GRAPH_GRAPHIO_HPP
 
