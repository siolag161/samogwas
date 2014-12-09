/****************************************************************************************
 * File: FLTM.cpp
 * Description: Implementation of a sequential version of the FLTM algorithm.
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 03/04/2014

 ***************************************************************************************/

#include "fltm/fltm.hpp"
#include "fltm/graph.hpp"
#include "fltm/latent_var_criteria.hpp"
#include "utils/matrix_utils.hpp"

using namespace utility;

namespace samogwas {

/**
 *
 */
void FLTM::operator()( FLTM_Result &result,
                       FLTM_Data &input,
                       FLTM_Options &options) {
        
  Matrix2GraphIndex mat2GraphIndex; StrLabel2GraphIndex label2GraphIndex;
  const unsigned nbrObsVars = input.matrix->size();
  std::cout << "\n========================================  beginning FLTM  ====================================" << std::endl;
  
  Graph &graph = *result.graph;
  initializeFLTM( input, result, mat2GraphIndex, label2GraphIndex, nbrObsVars, options.cardinality );
  
  for ( int step = 0; step < options.nbrSteps; ++step) {
    printf( "\n++++++++++++++++++++++++++++++++++++++++++++++++ beging step: %d of %d  ++++++++++++++++++++++++++++++++++++++++++++++++\n", step, options.nbrSteps );
    MatrixPtr nextStepMatrix(new Matrix); Matrix2GraphIndex nextStepMat2GraphIndex;
    if ( step > 0 )
      clustAlgo->invalidate();
    // clustAlgo->setData(input.matrix);
    printf( "start clustering of data of: %zu variables...\n", input.matrix->size() );
    auto partition = clustAlgo->run();

    // printf( "vares...\n", input.matrix->size() );

    auto clustering = partition.to_clustering(); // auto: to prevent explicite declaration of the type.
    printf( "end clustering, obtained: %zu clusters...\n", clustering.size() );

    // if (step==1) {
    //   for (auto &clt: clustering) {
    //     printf("cluster: ");
    //     for ( auto cit: clt ) {
    //       printf("%d ", cit);
    //     }
    //     std::cout << std::endl;
    //   }
    //   // exit(-1);
    // }
    int singletonCount = 0;
    if ( containsOnlySingletons( singletonCount, clustering) ) {
      std::cout << "stops due to only singleton. " << std::endl;
      return;
    }
     
    int goodClusterCount = 0; nextStepMatrix->reserve(clustering.size());
    int num_Clust = 0;
    for ( auto &cluster: clustering ) {
      std::cout << "processing cluster of sz: " << cluster.size() << std::endl;           

      if ( cluster.size() > 1 ) {
        MatrixPtr emMat(new Matrix); ResultEM resultEM;
        Variables emVars;
        std::cout << "creating latent var" << std::endl;           

        Variable latentVar = createLatentVar( boost::num_vertices(graph), cardFunc(cluster) );
        std::cout << "initing em" << std::endl;       

        initializeEM( *emMat, emVars, input, graph, cluster, mat2GraphIndex);
        std::cout << "Cluster: " << num_Clust++ << " over " << clustering.size() - singletonCount << " of step: " << step << std::endl;
        std::cout << "executing EM on setLabel of nbrVariables: " << cluster.size() << std::endl;
        emFunc->run( resultEM, latentVar, emVars, emMat, options.emThres );
        std::cout << "done EM" << std::endl << std::endl;        
        if ( goodLatentVariable( resultEM.imputedData, *input.matrix, cluster, options.latentVarQualityThres) ) {
          // std::cout << "it's a good cluster" << std::endl;    
          vertex_t vertex = addLatentNode( graph, latentVar, cluster, resultEM, label2GraphIndex ); // vertex is the index of the newly added node
          goodClusterCount++;
          nextStepMatrix->push_back(resultEM.imputedData);
          nextStepMat2GraphIndex.push_back( (int)vertex );
          result.imputedData->push_back(resultEM.imputedData );
          // std::cout << "moves on" << std::endl;    
        } else {
          // std::cout << "it's a bad cluster. init..." << std::endl;   
          initializeNextStep(*nextStepMatrix, nextStepMat2GraphIndex, mat2GraphIndex, *input.matrix, cluster);
          // std::cout << "i know, it's a bad cluster. what's next..." << std::endl << std::endl;   

        }
      } else {
        // std::cout << "it's a singleton luster. init..." << std::endl;           
        initializeNextStep(*nextStepMatrix, nextStepMat2GraphIndex, mat2GraphIndex, *input.matrix, cluster);
        // std::cout << "i know, it's a singleton cluster. what's next..." << std::endl << std::endl;   
      }
    }

    std::cout << "alrighty." << std::endl;           

    if ( goodClusterCount == 0 ) {
      std::cout << "stops due to 0 good clusters. " << std::endl;
      return;
    }
    mat2GraphIndex = nextStepMat2GraphIndex;
    *input.matrix = *nextStepMatrix;
    input.positions = extractPositionsForMatrixVariables( graph, mat2GraphIndex );
    printf( "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++ END  ++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", step, options.nbrSteps );

  }
}


///////////////////////////////////////////
void FLTM::initializeFLTM( FLTM_Data &input,
                           FLTM_Result &result,
                           Matrix2GraphIndex &mat2GraphIndex,
                           StrLabel2GraphIndex &label2GraphIndex,
                           const size_t &nbrVars,
                           const size_t &cardinality) {
  // `labels` stores the names (labels) of the variables.
  assert( input.matrix->size() == input.positions.size()); // `input.positions` has to be set.
  for ( size_t i = 0; i < input.labels.size(); ++i) {
    int level = 0;
    vertex_t vertex = createVertex( *result.graph, cardinality, true, input.labels[i], input.positions[i], level );
    label2GraphIndex[ input.labels[i] ] = (int) vertex;
    mat2GraphIndex.push_back((int)vertex) ;
  }
}

///////////////////////////////////////////////
std::vector<Position> FLTM::extractPositionsForMatrixVariables( const Graph &graph, Matrix2GraphIndex &mat2GraphIndex ) {
  std::vector<Position> localPositions( mat2GraphIndex.size(), 0 );
  for ( size_t i = 0; i < mat2GraphIndex.size(); ++i ) {
    localPositions[i] = graph[ mat2GraphIndex.at(i) ].position;
  }
  return localPositions;
}

///////////////////////////////////////
bool FLTM::containsOnlySingletons( int &singletonCount, const Clustering &clustering ) {
  assert(singletonCount == 0);
  for ( auto clt: clustering ) {      
    if ( clt.size() <= 1) {
      ++singletonCount;
    }
  }
  printf("there are: %d singleton clusters\n", singletonCount);
  return (clustering.size() == singletonCount);
}

///////////////////////////////////////
Variable FLTM::createLatentVar( const int lab, const int cardinality ) {
  return Variable( boost::lexical_cast<std::string>(lab), plIntegerType(0, cardinality - 1) );
}

///////////////////////////////////////
void FLTM::initializeEM( Matrix &emMat,
                         Variables &vars,
                         const FLTM_Data &input,
                         const Graph &graph,
                         const std::vector<int> &cluster,
                         const std::vector<int> &mat2GraphIndex )  {
  MatrixPtr tEMMat(new Matrix);
  tEMMat->reserve(cluster.size());
  vars.clear();
    // A column in `matrix` corresponds to an observation (e.g. an individual).
  // printf("start initing EM - 1 \n");
  tEMMat->push_back( std::vector<int>( ncols(*input.matrix), -1) ); // data for the latent variable (initialized to -1)
  // printf("start initing EM - 2 \n");
  for ( auto it: cluster ) {
    // printf("initing EM - %d vs %d\n", it, input.matrix->size());

    vars ^= graph[mat2GraphIndex.at(it)].variable;
    tEMMat->push_back( input.matrix->at(it) );      
  }
  // printf("start initing EM - 3 \n");

  emMat = Transpose(*tEMMat); // A row in `emMat corresponds to an observation (e.g. an individual).
  // delete(tEMMat);
  // printf("end initing EM\n\n");

}
///////////////////////////////////////////////
bool FLTM::goodLatentVariable( std::vector<int> &latentData,
                               Matrix &emMat,
                               std::vector<int> &cluster,
                               double latentVarQualityThres )
{
  AverageMutInfo averageMutInfo;
  double measuredQuantity = averageMutInfo( latentData, emMat, cluster );
  // printf("evaluting quality: %f vs threshold: %f\n", measuredQuantity, latentVarQualityThres);
  return ( measuredQuantity >= latentVarQualityThres);
}

///////////////////////////////////////////////
vertex_t FLTM::addLatentNode( Graph &graph,
                              const Variable &latentVar,
                              std::vector<int> &cluster,
                              ResultEM &resultEM,
                              StrLabel2GraphIndex &label2GraphIndex ) {
  const vertex_t vertex = createVertex( graph,
                                        latentVar.cardinality(),
                                        false, // isLeaf = false
                                        latentVar.name() );
  graph[vertex].setupProperties(&graph, resultEM.jointDistribution, label2GraphIndex);
  for ( auto clt: cluster ) {
    boost::add_edge(vertex, clt, graph);
  }
  label2GraphIndex[latentVar.name()] = (int) vertex;
  return vertex;
}

/////////////////////////////////////////////////////////////////////////
/** nextStepMatrix: The newly created latent variables have been added.
  *                 All the children of these latent variables have been discarded.
  *
  */
void FLTM::initializeNextStep( Matrix &nextStepMatrix,
                               Matrix2GraphIndex &nextStepMat2GraphIndex,
                               const Matrix2GraphIndex &mat2GraphIndex,
                               const Matrix &matrix,
                               const std::vector<int> &cluster) {
  // printf("alibaba-1\n");
  for ( auto cit: cluster ) {
    nextStepMat2GraphIndex.push_back( mat2GraphIndex[cit] );
    nextStepMatrix.push_back( matrix[cit] );
  }
  // printf("alibaba-3-done\n");

}

} // namespace samogwas ends here.

 
