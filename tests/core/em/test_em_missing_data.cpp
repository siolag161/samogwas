#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE
#endif  
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <map>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <math.h>
#include <boost/lexical_cast.hpp>
#include "em/em.hpp"

#include "data_generation.hpp"
#include "clustering/cast.hpp"
#include "clustering/dbscan.hpp"
#include "distance/dissimilarity.hpp"
#include "fltm/graph.hpp"
#include "utils/matrix_utils.hpp"
#include "em/core_em.hpp"
#include <initializer_list>

#define NOT_MISSING 1
#define MISSING 0
#define MISSING_VALUE -1
class Data 
{ 
};

typedef samogwas::MutInfoDissimilarity<data_gen::Matrix> MutInfoDiss;
typedef samogwas::MutInfoSimilarity<data_gen::Matrix> MutInfoSimi;
typedef samogwas::DBSCAN<MutInfoDiss> DBSCAN;
typedef samogwas::CAST<MutInfoSimi> CAST;
typedef samogwas::Partition Partition;

using namespace data_gen;// using namespace samogwas;
using namespace utility;
typedef std::vector< std::vector<bool>> DefTab;
typedef std::shared_ptr<DefTab> DefTabPtr;
template<typename T>
void printMatrix(T& mat);


BOOST_FIXTURE_TEST_SUITE( Test_EM_Missing_Data, Data ) 

inline void prepareEM( MatrixPtr& emMat,
                       Variables& vars,
                       const MatrixPtr cltMat,
                       const samogwas::Graph& graph,
                       const std::vector<int>& cluster,
                       const std::vector<int> local2Global);

inline DefTabPtr getDefTab(Matrix& mat);

Variable createLatentVar( const int lab, const int cardinality ) {
  return Variable( boost::lexical_cast<std::string>(lab), plIntegerType(0, cardinality - 1) );
}

BOOST_AUTO_TEST_CASE( Test_EM_functional ) {
  
  size_t nclusts = 5, ncols = 40;
  size_t N = 3, CARD = 3, MAX_POS = 50;
  int nrows = nclusts*N;
  std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  auto data = GenerateClusteredData( nclusts, N, CARD, ncols )();  
  for ( int i = 0; i < 15; ++i ) {
    (*data)[i][i] = MISSING_VALUE;
  }

  // std::cout << "dim of mat: "
  // printf("dim of data: %d-%d\n", (int)data->size(), (int)(*data)[0].size());
  
  auto simi = std::make_shared<MutInfoSimi>(data, positions, MAX_POS, MISSING_VALUE);
  
  CAST cast( simi, 0.3 );  
  Partition result = cast();
  
  for ( int i = 0; i < nrows; ++i ) {
    int expected_cluster = i / 3;
    BOOST_CHECK_EQUAL(result.getLabel(i), expected_cluster );
  }

  std::vector<int> local2Global;
  samogwas::Graph graph;
  for (size_t i = 0; i < N*nclusts; ++i) {
    samogwas::vertex_t vertex = createVertex( graph, 3, true,
                                    boost::lexical_cast<std::string>(i), i, 0 );

    local2Global.push_back(vertex) ;
  }

  int i = 0;
  
  for ( auto& clt: result.to_clustering() ) {
    auto emMat = std::make_shared<Matrix>();
    Variables vars;
    prepareEM( emMat, vars, data, graph, clt, local2Global );
    Variable latentVar = createLatentVar( boost::num_vertices(graph), 3 );
    samogwas::ResultEM resultEM;
    samogwas::NaiveBayesEM multiEM(3,1);
    auto defTab = samogwas::EMInterface::createDefinitionTable(emMat);
    // printf("some dim (%ul, %ul) vs (%ul, %ul)\n",
    //        utility::nrows(*emMat), 
    //        utility::ncols(*emMat),
    //        utility::nrows(*defTab), 
    //        utility::ncols(*defTab) );
    // multiEM( resultEM, latentVar, vars, emMat, 0.000001, defTab );

    // printMatrix(emMat);
    // printMatrix(*defTab);

    // for (int i=0; i<emMat.size();++i) {
    //   for (int j=0; j<emMat[0].size();++j) {
    //     // std::cout << "("<<(*emMat)[i][j]<< " - " << (*defTab)[i][j] << "), "; // 
    //   }
    //   // std::cout<<std::endl;
    // }
    
    // break;
  }
}

BOOST_AUTO_TEST_CASE( Test_EM_bool ) {

   auto emMat = std::make_shared<Matrix>(std::initializer_list<std::vector<int>>{
      std::vector<int>{ MISSING_VALUE,1,1,1},
      std::vector<int>{ MISSING_VALUE,0,1,2},
      std::vector<int>{ MISSING_VALUE,0,1,2},
      std::vector<int>{ MISSING_VALUE,2,2,0} });
  auto defTab = std::make_shared<DefTab>(std::initializer_list< std::vector<bool> >{
    std::vector<bool>{MISSING, NOT_MISSING, NOT_MISSING, NOT_MISSING},
    std::vector<bool>{MISSING, NOT_MISSING, NOT_MISSING, NOT_MISSING},
    std::vector<bool>{MISSING, MISSING, NOT_MISSING, NOT_MISSING},
    std::vector<bool>{MISSING, NOT_MISSING, NOT_MISSING, NOT_MISSING}    
    });


  samogwas::ResultEM resultEM;
  samogwas::NaiveBayesEM multiEM(3,1);
  Variable latentVar = Variable( "3", plIntegerType(0, 2) );
  Variables vars;
  vars ^= Variable( "0", plIntegerType(0, 2) );
  vars ^= Variable( "1", plIntegerType(0, 2) );
  vars ^= Variable( "2", plIntegerType(0, 2) );

  // multiEM( resultEM, latentVar, vars, emMat, 0.000001, defTab );
  
}


inline void prepareEM( MatrixPtr& emMat,
                Variables& vars,
                const MatrixPtr cltMat,
                const samogwas::Graph& graph,
                const std::vector<int>& cluster,
                const std::vector<int> local2Global ) {
  
  auto tEMMat = std::make_shared<Matrix>();
  tEMMat->reserve(cluster.size());
  vars.clear();  
  tEMMat->push_back( std::vector<int>( ncols(*cltMat), MISSING_VALUE) );  
  for ( auto& it: cluster ) {    
    vars ^= graph[local2Global.at(it)].variable;
    tEMMat->push_back( cltMat->at(it) );      
  }
  emMat = Transpose(*tEMMat);
   // delete(tEMMat);
}


DefTabPtr getDefTab(MatrixPtr dataMat) {
  const auto nbrInds = utility::nrows(*dataMat);
  const auto nbrVars = utility::ncols(*dataMat);
  // printf("size: %d,%d\n", nbrInds, nbrVars);
  auto defTable =
      std::make_shared<DefTab>( nbrInds, std::vector<bool>(nbrVars, NOT_MISSING)) ;
  for (size_t ind = 0; ind < nbrInds; ++ind) {
    // defTable.push_back(std::vector<bool>(nbrVars, NOT_MISSING));
    for (size_t var = 1; var < nbrVars; ++var) {
      if ((*dataMat)[ind][var] == DATA_MISSING_VALUE) {
        //   // printf("(%d,%d): ", ind, var);                
        //   // for (size_t i = 0; i < nbrInds; ++i) {
        (*defTable)[ind][var] = 0;
        //   // }
        //   // printf("(%d,%d): ", ind, var);
        //   // std::cout << defTable[ind][var] << std::endl;
        //   // printf
      }
    }
    // std::cout << std::endl;
    (*defTable)[ind][0] = MISSING; // The first column is FALSE.
  }

  return defTable;  
}


// template<typename M>
// void printMatrix(M& mat) {
//   for (int i=0; i<mat.size();++i) {
//     for (int j=0; j<mat[0].size();++j) {
//       std::cout<<mat[i][j]<< ", ";
//     }
//     std::cout<<std::endl;
//   }
// }


BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here
