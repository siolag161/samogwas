/****************************************************************************************
 * File: community_generation.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 21/11/2014

 ***************************************************************************************/
#ifndef TEST_COMMUNITY_GENERATION_HPP
#define TEST_COMMUNITY_GENERATION_HPP

#include "data_generation.hpp"
#include "clustering/louvain/graph.hpp"

samogwas::louvain::Graph* generateGraph(int commCard, int commCount);
static size_t N = 3, CARD = 3, MAX_POS = 5, NCOLS = 40;

#include "distance/dissimilarity.hpp"
#include "distance/similarity.hpp"

typedef samogwas::MutInfoDissimilarity<data_gen::Matrix> MutInfoDiss;
typedef samogwas::MutInfoSimilarity<data_gen::Matrix> MutInfoSimi;

/****************************************************************************************/
#endif // _COMMUNITY_GENERATION_HPP
