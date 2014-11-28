#include "data_generation.hpp"
#include "community_generation.hpp"

using namespace samogwas;
using namespace samogwas::louvain;
using namespace data_gen;

#include "distance/dissimilarity.hpp"
#include "distance/similarity.hpp"


Graph* generateGraph(int commCard, int commCount) {
  // Graph g;

  // size_t nclusts = 5, ncols = 40;
  //size_t N = 3, CARD = 3, MAX_POS = 50;
  int nrows = commCard*commCount;
  std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  auto data = GenerateClusteredData( commCount, commCard, CARD, NCOLS )();  
  //  std::shared_ptr<MutInfoSimi> simi(new MutInfoSimi(data, positions, MAX_POS, -1));
  // SimilarityMatrix* simi = new MutInfoSimi(data, positions, MAX_POS, -1);
  std::shared_ptr<SimilarityMatrix> simi(new MutInfoSimi(data, positions, MAX_POS, -1));

  Graph* g  =  new Graph(simi);
  // for (int i = 0; i<3;++i) {
  //   for (int j = 0; j<3;++j) {
  //     std::cout << g->weight(i,j) << std::endl;
  //   }
  // }
  std::cout << "---------------" << std::endl;

  // Graph* g2  =  g;
  // for (int i = 0; i<10;++i) {
  //   for (int j = 0; j<10;++j) {
  //     std::cout << g2->weight(i,j) << std::endl;
  //   }
  // }
  std::cout << "---------------" << std::endl;

  // Graph* g2 = g;
  // for (int i = 0; i<10;++i) {
  //   for (int j = 0; j<10;++j) {
  //     std::cout << g2->weight(i,j) << std::endl;
  //   }
  // }
  // std::cout << g.weight(1,5) << std::endl;
  
  return g;

  // return g;
}
