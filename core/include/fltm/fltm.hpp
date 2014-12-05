/****************************************************************************************
 * File: FLTM.hpp
 * Description: This module provides the implementation for the construction of the FLTM (Forest of Latent Tree Models).
 * -----------: This module requires an input data matrix to be row-major ( a row represents a variable, a column
 * -----------: represents an observation (e.g. an individual in the case of GWASs)).
 * @ref: @todo: BMC Bioinformatics
 * 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 09/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_FLTM_HPP
#define SAMOGWAS_FLTM_HPP  

#include <vector>
#include <memory>

#include "clustering/clustering.hpp"
#include "em/em.hpp"
#include "utils/matrix_utils.hpp"

#include "em_card_func.hpp"
#include "graph.hpp"
#include "graph_io.hpp"
#include "core_fltm.hpp"

namespace samogwas
{
typedef std::shared_ptr<EMInterface> EM_Algo_Ptr;
typedef std::shared_ptr<AlgoClusteringInterface> Clust_Algo_Ptr;
struct FLTM {

  FLTM( Clust_Algo_Ptr clustA,
        CardFunc &cardF,
        EM_Algo_Ptr emF): clustAlgo(clustA), cardFunc(cardF), emFunc(emF) { }

  void operator()( FLTM_Result &result, FLTM_Data &data, FLTM_Options &options);
    
  // ~FLTM() { delete clustAlgo; delete emFunc; }

 protected:
  void initializeFLTM( FLTM_Data &input,
                       FLTM_Result &result,
                       Matrix2GraphIndex &mat2GraphIndex,
                       StrLabel2GraphIndex &label2GraphIndex,
                       const size_t &nbrVars,
                       const size_t &cardinality);

  std::vector<Position> extractPositionsForMatrixVariables( const Graph &graph, Matrix2GraphIndex &mat2GraphIndex  );

  bool containsOnlySingletons( int &singletonCount,
                               const Clustering &clustering );

  Variable createLatentVar( const int lab, const int cardinality );

  ///////////////////////////////////////
  void initializeEM( Matrix &emMat,
                     Variables &vars,
                     const FLTM_Data &input,
                     const Graph &graph,
                     const std::vector<int> &cluster,
                     const Matrix2GraphIndex &mat2GraphIndex );

  bool goodLatentVariable( std::vector<int> &latentData,
                           Matrix &globalMatrix,
                           std::vector<int> &cluster,
                           double latentVarQualityThres );

  vertex_t addLatentNode( Graph &graph,
                          const Variable &latentVar,
                          ResultEM &resultEM,
                          StrLabel2GraphIndex &label2GraphIndex );

  void initializeNextStep( Matrix &nextStepMatrix, Matrix2GraphIndex &nextStepMat2GraphIndex,
                           const Matrix2GraphIndex &mat2GraphIndex,
                           const Matrix &globalMatrix, const std::vector<int> &cluster);

//  void initializeGraph( const FLTM_Data &input, Graph &graph );

 protected:
  Clust_Algo_Ptr clustAlgo;
  CardFunc &cardFunc;
  EM_Algo_Ptr emFunc;
};

} // namespace samogwas ends here.



/****************************************************************************************/
#endif //SAMOGWAS_FLTM_HPP
