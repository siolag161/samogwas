/****************************************************************************************
 * File: core_fltm.hpp
 * Description: All FLTM algorithms must implement this interface.
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 01/08/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_CORE_FLTM_HPP
#define SAMOGWAS_CORE_FLTM_HPP

#include <vector>
#include <map>

#include "clustering/clustering.hpp"
#include "em/em.hpp"
#include "utils/matrix_utils.hpp"

#include "em_card_func.hpp"
#include "graph.hpp"


namespace samogwas
{

typedef std::vector< std::vector<int> > Matrix;
typedef plSymbol Variable;
typedef int Index;
typedef std::string StrLabel;
typedef int Position;
typedef std::map<StrLabel, Index> StrLabel2GraphIndex; // for a variable, label to its node index in the graph
// typedef std::map<StrLabel, Position> StrLabel2Pos;
typedef std::vector<Index> Matrix2GraphIndex; // maps a variable (i.e. row) in a matrix to its node index in the graph.

/**
 *
 */
struct FLTM_Result {
  // FLTM_Result(): nbrLatentVariables(0) {}

  // void addNode(const Node& node) {
  //   while( (node.level) >= level2LatentVars.size()) { // eventually creates empty levels.
  //     level2LatentVars.push_back(std::vector<vertex_t>());
  //   }
  //   level2LatentVars[node.level].push_back(node.index);
  //   ++nbrLatentVariables;

  //   Graph graph;
  // }

  // ///
  // int nbrLatentVariables;
  // std::vector< std::vector<vertex_t> > level2LatentVars; // latent variables by level
  Matrix imputedData;
  Graph graph;    
};

/** Encapsulates the information and the structure (graph) of the input data.
 *
 */
struct FLTM_Data {
  std::vector<StrLabel> labels;
  std::vector<Position> positions; // physical positions
  std::vector<unsigned> indexes; //
  std::vector< std::vector<int> > matrix;
  // Graph graph;
  int cardinality;
};

struct FLTM_Options {
  int cardinality;
  int nbrSteps;
  double emThres; // controls EM algorithm convergence.
  double latentVarQualityThres;
};


} // namespace samogwas ends here.

/****************************************************************************************/
#endif // SAMOGWAS_CORE_FLTM_HPP
