/****************************************************************************************
 * File: core_fltm.hpp
 * Description: 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 24/08/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_CORE_FLTM_HPP
#define SAMOGWAS_CORE_FLTM_HPP

namespace samogwas
{


typedef std::vector< std::vector<int> > Matrix;
typedef plSymbol Variable;
typedef int Index;  
typedef std::string StrLabel;
typedef int Position;
typedef std::map<StrLabel, Index> StrLabel2GraphIndex; // for a variable, label to global index
typedef std::map<StrLabel, Position> StrLabel2Pos;
typedef std::vector<Index> Matrix2GraphIndex;

struct FLTM_Result {
  FLTM_Result(): nbrLatentVariables(0) {}  
  void addNode(const Node& node) {
    while( (node.level) >= level2LatentVars.size()) {
      level2LatentVars.push_back(std::vector<vertex_t>());
    }
    level2LatentVars[node.level].push_back(node.index);
    ++nbrLatentVariables;
  }
  int nbrLatentVariables;
  std::vector< std::vector<vertex_t> > level2LatentVars; 
  Matrix imputedData;
};


struct FLTM_Data {  
  std::vector<StrLabel> labels;
  std::vector<Position> positions;
  std::vector<unsigned> ids;
  std::vector< std::vector<int> > matrix;
  Graph graph;
  int cardinality;
};

struct FLTM_Options {
  int cardinality;
  int nbrSteps;
  double emThres;
  double infoThres;
};


} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_CORE_FLTM_HPP
