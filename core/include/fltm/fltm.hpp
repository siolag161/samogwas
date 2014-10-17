/****************************************************************************************
 * File: FLTM.hpp
 * Description: 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 09/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_FLTM_HPP
#define SAMOGWAS_FLTM_HPP  

#include <vector>
#include "clustering/clustering.hpp"
#include "em/em.hpp"
#include "utils/matrix_utils.hpp"

#include "em_card_func.hpp"
#include "graph.hpp"
#include "graph_io.hpp"
#include "core_fltm.hpp"

namespace samogwas
{
 
struct FLTM {
  FLTM( AlgoClusteringInterface* clustA,
        CardFunc& cardF,
        EMInterface* emF): clustAlgo(clustA), cardFunc(cardF), emFunc(emF) { }

  void operator()( FLTM_Result& result, FLTM_Data& data, FLTM_Options& opt );
  ~FLTM() { delete clustAlgo; delete emFunc; }

 protected:
  void initializeFLTM( FLTM_Data &input,
                       FLTM_Result &result,
                       Matrix2GraphIndex &mat2GraphIndex,
                       StrLabel2GraphIndex &label2GraphIndex,
                       const size_t &nbrVars,
                       const size_t &cardinality);

  std::vector<Position> extractPositionsForMatrixVariables( const Graph& graph, Matrix2GraphIndex& mat2GraphIndex );

  bool containsOnlySingletons( int& singleton,
                               const Clustering& clustering );

  Variable createLatentVar( const int lab, const int cardinality );

  ///////////////////////////////////////
  void initializeEM( Matrix &emMat,
                     Variables &vars,
                     const FLTM_Data &input,
                     const Graph& graph,
                     const std::vector<int> &cluster,
                     const std::vector<int> local2Global);

  bool goodLatentVariable( std::vector<int>& latentCol,
                           Matrix& transposedMat,
                           std::vector<int>& cluster,
                           double goodLatentVarThres );

  vertex_t addLatentNode( Graph& graph,
                          const Variable& latentVar,
                          ResultEM& resultEM,
                          StrLabel2GraphIndex& label2GraphIndex );

  void initializeNextStep( Matrix &nextRowMatrix, Matrix2GraphIndex &nextRoundMat2GraphIndex,
                           const Matrix2GraphIndex &mat2GraphIndex,
                           const Matrix &matrix, const std::vector<int> &cluster);

  void initializeGraph( const FLTM_Data &input, Graph& graph );

 protected:
  AlgoClusteringInterface* clustAlgo;
  CardFunc& cardFunc;
  EMInterface* emFunc;
};

} // namespace samogaws ends here.



/****************************************************************************************/
#endif // end SAMOGWAS_FLTM_HPP
