/****************************************************************************************
 * File: clustering_parameters.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 05/12/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_COMMON_CLUSTERING_PARAMETERS_HPP
#define SAMOGWAS_COMMON_CLUSTERING_PARAMETERS_HPP

#include "clustering/clustering.hpp"
#include "distance/dissimilarity.hpp"
#include "distance/similarity.hpp"
#include "distance/comparable.hpp"

#include "clustering/clustering.hpp"
#include "clustering/cast.hpp"
#include "clustering/dbscan.hpp"
#include "clustering/louvain/louv.hpp"

#include <vector>
#include <memory>
#include <stdexcept>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>

namespace samogwas
{

typedef AlgoClusteringInterface Algorithm;
typedef std::shared_ptr<Algorithm> ClustAlgoPtr;
typedef std::vector< std::vector<int> > Matrix; // We consider here only vector of int is relevant
typedef int Position;
typedef MutInfoSimilarity<Matrix> MutInfoSimi;
typedef MutInfoDissimilarity<Matrix> MutInfoDiss;
typedef std::vector<Position> Positions;
typedef std::shared_ptr<Matrix> MatrixPtr;

typedef DBSCAN<MutInfoDiss> DBSCAN_Algo;
typedef CAST<MutInfoSimi> CAST_Algo;

// inline AlgoClusteringInterface* getAlgoClust( FLTM_Data& input, Options& opt );

inline std::vector<ClustAlgoPtr> read_clustering_algos( MatrixPtr matrix, Positions& positions,
                                                        double maxDist, double simiThres, std::istream& is) {
  using boost::property_tree::ptree;
  ptree pt;
  read_xml(is, pt);

  std::vector<ClustAlgoPtr> rs;
  for( const ptree::value_type &v: pt.get_child("clustering") ) {
    if( v.first == "algorithm" ) {
      ClustAlgoPtr algorithm;

      auto algo_cfg = v.second;
      if ( algo_cfg.get<std::string>("name") == "DBSCAN" ) {
        auto diss = std::make_shared<MutInfoDiss>( matrix, positions, maxDist, simiThres );

        
        int minPts; double eps;
        for( const ptree::value_type &pam: algo_cfg.get_child("parameters") ) {
          auto dat = pam.second.data();
          if ( pam.second.get<std::string>("<xmlattr>.name") == "minPts") {
            minPts = boost::lexical_cast<int>(dat);
          } else if (pam.second.get<std::string>("<xmlattr>.name") == "eps") {
            eps = boost::lexical_cast<double>(dat);
          }
        }
        algorithm = std::make_shared<DBSCAN_Algo>( diss, minPts, eps );      
      } else if ( algo_cfg.get<std::string>("name") == "CAST" ) {
        auto simi = std::make_shared<MutInfoSimi>( matrix, positions, maxDist, simiThres );
        double cast;
        for( const ptree::value_type &pam: algo_cfg.get_child("parameters") ) {
          auto dat = pam.second.data();
          if ( pam.second.get<std::string>("<xmlattr>.name") == "t" ) {
            cast = boost::lexical_cast<double>(dat);
          }
        }
        algorithm = std::make_shared<CAST_Algo>( simi, cast);

      } else if ( algo_cfg.get<std::string>("name") == "LOUVAIN" ) {
        auto simi = std::make_shared<MutInfoSimi>( matrix, positions, maxDist, simiThres );           
        algorithm = std::make_shared<louvain::MethodLouvain>(simi);
      } else {
        throw std::invalid_argument( "unknown algorithm" );
      }
      printf("processing: %s\n", algorithm->name());
      rs.push_back(algorithm);
    }
  }    
  return rs;
}

} // namespace samogwas ends here. 

/****************************************************************************************/
#endif // SAMOGWAS_CLUSTERING_PARAMETERS_HPP
