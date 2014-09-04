/****************************************************************************************
 * File: FLTM.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 09/07/2014

 ***************************************************************************************/
#ifndef FLTM_SAMOGWAS_HPP
#define FLTM_SAMOGWAS_HPP  

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
  FLTM( AlgoClustering* clustA, CardFunc& cardF, EMFunc* emF): clustAlgo(clustA), cardFunc(cardF), emFunc(emF) { }  

  void operator()( FLTM_Result& result, FLTM_Data& data, FLTM_Options& opt );
  ~FLTM() { delete clustAlgo; delete emFunc; }

 protected:
  void setupVariables( FLTM_Data& input,
                       FLTM_Result& result,
                       Matrix2GraphIndex& mat2GraphIndex,
                       Label2GraphIndex& label2GraphIndex,
                       const size_t& nbrVars,
                       const size_t& cardinality);

  std::vector<Position> getLocalPositions( const Graph& graph, Matrix2GraphIndex& mat2GraphIndex  );

  bool containsOnlySingletons( int& singleton,
                               const Clustering& clustering );
  Variable createLatentVar( const int lab, const int cardinality );
  ///////////////////////////////////////
  void prepareEM( Matrix& emMat,
                  Variables& vars,
                  const FLTM_Data& input,
                  const std::vector<int>& cluster,
                  const std::vector<int> local2Global );

  bool goodLatentVariable( std::vector<int>& latentCol,
                           Matrix& transposedMat,
                           std::vector<int>& cluster,
                           double goodLatentVarThres );

  vertex_t addLatentNode( Graph& graph,
                          const Variable& latentVar,
                          ResultEM& resultEM,
                          Label2GraphIndex& label2GraphIndex );

  void updateNextRow( Matrix& nextRowMatrix, Matrix2GraphIndex& nextRoundMat2GraphIndex,
                      const Matrix2GraphIndex& mat2GraphIndex, const Matrix& matrix, const std::vector<int>& cluster );

 protected:
  AlgoClustering* clustAlgo;
  CardFunc& cardFunc;
  EMFunc* emFunc;
};

} // namespace fltmends here. fltm



/****************************************************************************************/
#endif // FLTM_FLTM_HPP
