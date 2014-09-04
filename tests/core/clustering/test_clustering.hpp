/****************************************************************************************
 * File: test_clustering.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 30/07/2014

 ***************************************************************************************/
#ifndef _TEST_CLUSTERING_HPP
#define _TEST_CLUSTERING_HPP

#include "data_generation.hpp"
#include "clustering/cast.hpp"
#include "clustering/dbscan.hpp"
#include "distance/dissimilarity.hpp"

typedef samogwas::MutInfoDissimilarity<data_gen::Matrix> MutInfoDiss;
typedef samogwas::MutInfoSimilarity<data_gen::Matrix> MutInfoSimi;

typedef samogwas::DBSCAN<MutInfoDiss> DBSCAN;
typedef samogwas::CAST<MutInfoSimi> CAST;  


typedef samogwas::Partition Partition;

using namespace data_gen;

/****************************************************************************************/
#endif // _TEST_CLUSTERING_HPP
