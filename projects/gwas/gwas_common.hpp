/****************************************************************************************
 * File: gwas_common.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 10/12/2014

 ***************************************************************************************/
#ifndef _GWAS_COMMON_HPP
#define _GWAS_COMMON_HPP

#include "fltm/graph.hpp"
#include "fltm/graph_io.hpp"
#include "app_options.hpp"

////////////////////////////////////////////////////////////////////
typedef std::vector<int> PhenoVec;
typedef std::shared_ptr<PhenoVec> PhenoVecPtr;
typedef std::vector<std::vector<double>> ValueMat;
typedef std::shared_ptr<ValueMat> ValueMatPtr;

ValueMatPtr loadScores(std::string& infile);
ValueMatPtr loadThres(std::string& infile);
PhenoVecPtr loadPhenotype(std::string& phenoFile);

////////////////////////////////////////////////////////////////////
void assureGraphPositions( samogwas::Graph& g );

char* current_date();

boost::filesystem::path outputDir( std::string& outputDir, bool date = true);

PhenoVecPtr loadPhenotype(std::string& phenoFile);

void count_cluster(Graph&g, std::vector<double>& scores);

std::vector<int> getGraphParent( const Graph& graph );

std::vector<int> countClusterSiblings( Graph& graph );

/****************************************************************************************/
#endif // _GWAS_COMMON_HPP
