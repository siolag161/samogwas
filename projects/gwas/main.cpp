// /****************************************************************************************
//  * File: fltm.cpp
//  * Description: 
//  * @author: siolag161 (thanh.phan@outlook.com)
//  * @date: 09/07/2014
//  ***************************************************************************************/

#include <iostream>
#include <thread>
#include <chrono>
#include <memory>
#include <limits>

#include <fstream>
#include <cstdio>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp> // to obtain the program's name

#include "distance/comparable.hpp"

#include "clustering/cast.hpp"
#include "clustering/dbscan.hpp"
#include "clustering/louvain/louv.hpp"

#include "gwas/gwas_basic_strategy.hpp"

#include "distance/dissimilarity.hpp"
#include "distance/similarity.hpp"
#include "fltm/fltm.hpp"

#include "utils/option_printer.hpp"
#include "utils/custom_option_desc.hpp"

#include "utils/logs_utils.hpp"
#include "data_load.hpp"
#include "app_options.hpp"

#include "statistics/permutation_test.hpp"
#include "gwas/gwas_basic_strategy.hpp"
#include "gwas/score_node_criteria.hpp"

#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <memory>
#include <boost/lockfree/queue.hpp>

#include "gwas_common.hpp"

using namespace utility;
using namespace samogwas;

//////////////////////////////////////////////////////////////
Vertex addVirtualRoot( Graph& g, std::vector<Vertex>& parent);
void saveGWASSearching( std::string& outFileName, Graph& graph, std::vector<Vertex>& parent,
                 std::vector<Vertex>& visited, std::vector<std::vector<double>>& scores );
void saveRegions( std::string& outFileName, int chromosome, Graph& graph );

void assureGraphPositions(Graph& graph);

//////////////////////////////////////////////////////////////
int main( int argc, char** argv ) {
  auto pos = getToolProgramOptions( argc, argv );
  std::cout << "Loading graph data...\n" << std::endl;
  Graph graph;
  BayesGraphLoad graphLoad;
  graphLoad( graph,
             pos.inputLabelFile,
             pos.bayesVertices,
             pos.bayesDist );

  // assurePosition(graph);

  // printf("load graph of %d edges and %d vertices\n", boost::num_edges(*result.graph),
  //        boost::num_vertices(*result.graph));

  printf("Done loading graph of %d edges and %d vertices\n", boost::num_edges(graph),
         boost::num_vertices(graph));
  
  auto outputPath = outputDir(pos.outputDir);
  if ( pos.task == 0 ) {    
    FLTM_Data fltm_data;  
    std::cout << "Loading geno data from " << pos.inputDataFile << std::endl; // todo: logging
    fltm_data.matrix = loadDataTable(pos.inputDataFile);
    std::cout << "Loading pheno data from " << pos.inputPheno<< std::endl; // todo: logging
    auto pheno = loadPhenotype( pos.inputPheno );
    std::cout << "Loading label data from " << pos.inputLabelFile << std::endl; // todo: logging
    loadLabelPosition2( fltm_data.labels, fltm_data.indexes, fltm_data.positions, pos.inputLabelFile );
  
    auto scores = loadScores(pos.scoreFile);
    auto thres = loadThres(pos.thresFile);
    auto criteria = std::make_shared<MultiSourceNodeCriterion>(*thres, *scores, pos.overall_thres); 
  
    typedef std::map<Vertex, boost::default_color_type> ColorMap;
    typedef boost::associative_property_map<ColorMap> Color;
    typedef std::map<Vertex, double> ScoreMap;
    std::vector<Vertex> parent( boost::num_vertices(graph), -1);  
    auto root = addVirtualRoot(graph, parent);    
    ColorMap cmap;
    auto colorMap = std::make_shared<Color>(cmap);
    GWAS_Basic_Visitor visitor( criteria, std::make_shared<ScoreMap>(), colorMap);
    boost::queue<vertex_t> q;
    boost::breadth_first_search( graph, root, q, visitor, *colorMap);
    boost::remove_vertex(root,graph); 
    auto outFileName = ( outputPath / "gwas_result.txt").string();
    saveGWASSearching( outFileName, graph, parent, visitor.visitedVertices(), *scores );
  } else {
    auto outFileName = ( outputPath / "regions.txt").string();
    std::cout << "performing saving regions task to: " << outFileName << std::endl;
    saveRegions( outFileName, pos.chromosome, graph );   
  }
  std::cout << "done! program now exits\n";
}


ValueMatPtr loadScores(std::string& infile) {
  auto scores = std::make_shared<ValueMat>();

  std::ifstream scoreFile(infile.c_str());
  if (!scoreFile) {
    std::cout << "score file does not exist. program aborted\n";
    exit(-1);
  }
  
  utility::CSVIterator<double> scoreLine(scoreFile);
  for( ; scoreLine != utility::CSVIterator<double>(); ++scoreLine ) {
    std::vector<double> line(scoreLine->size(), 0.0);
    for (size_t i = 0; i < scoreLine->size(); ++i)
      line[i] = scoreLine->at(i);
    scores->push_back(line);
  }

  return scores;
}    
     
ValueMatPtr loadThres(std::string& infile) {
  auto thres = std::make_shared<ValueMat>();

  std::ifstream thresFile(infile.c_str());
  if (!thresFile) {
    std::cout << "threshold file does not exist. program aborted\n";
    exit(-1);
  }
  
  utility::CSVIterator<double> thresLine(thresFile);
  for( ; thresLine != utility::CSVIterator<double>(); ++thresLine ) {    
    std::vector<double> line(thresLine->size(), 0.0);
    for (size_t i = 0; i < thresLine->size(); ++i)
      line[i] = thresLine->at(i);
    thres->push_back(line);
  }

  return thres;
}

//////////////////////////////////////


Vertex addVirtualRoot( Graph& g, std::vector<Vertex>& parent ) {
  Vertex root = boost::add_vertex(g);

  boost::graph_traits<Graph>::vertex_iterator i, end;

  for ( auto vp = boost::edges(g); vp.first != vp.second; ++vp.first ) {
    auto s = boost::source(*vp.first,g);
    auto t = boost::target(*vp.first,g);
    parent[t] = s;
  }
    
  for ( auto vi = boost::vertices(g); vi.first != vi.second; ++vi.first ) {
    auto v = *vi.first;

    if ( parent[v] == -1 ) {
      boost::add_edge(root,v,g);
      printf("adding %d -> %d\n", root, v);
    }
    
  }

  return root;
}

////////////////////////////////////////

// static const std::string ID = "id"; static const std::string LABEL = "label";
// static const std::string LATENT = "latent"; static const std::string PARENT = "parent";
// static const std::string LEVEL = "level"; static const std::string POSITION = "position";
// static const std::string CARDINALITY = "cardinality";
// static const std::string PARENT_ID = "parent_id";
// static const char SEPARATOR = ',';


void saveGWASSearching( std::string& outFileName, Graph& graph, std::vector<Vertex>& parent,
                 std::vector<Vertex>& visited, std::vector<std::vector<double>>& scores ) {
  std::ofstream rsOut(outFileName);
  rsOut << ID << SEPARATOR
        << LATENT << SEPARATOR
        << PARENT << SEPARATOR
        << LEVEL << SEPARATOR
        << POSITION << SEPARATOR
        << CARDINALITY << "\n";  // writes header

  for ( auto v: visited ) {
    auto node = graph[v];
    auto latent = !node.isLeaf;
    rsOut << v << SEPARATOR
          << latent << SEPARATOR
          << graph[v].label << SEPARATOR
          << parent[v] << SEPARATOR
          << node.level << SEPARATOR
          << node.position << SEPARATOR
          << node.variable.cardinality() << std::endl;
          // << SEPARATOR
  }

  rsOut.close();

}

void saveRegions( std::string& outFileName, int chromosome, Graph& graph ) {
  const char SEPARATOR = ' ';
  std::ofstream rsOut(outFileName);
  for ( auto vi = boost::vertices(graph); vi.first != vi.second; ++vi.first ) {
    auto v = *vi.first;
    Node& node = graph[v];
    int sz_start = std::numeric_limits<int>::max(), sz_end = -std::numeric_limits<int>::max();

    int sum = 0, pos = 0, count = 0;
    for ( auto ei = boost::out_edges(v, graph); ei.first != ei.second; ++ei.first ) {
      auto child = graph[boost::target(*ei.first, graph)];
      sz_start = std::min(sz_start, child.position);
      sz_end = std::max(sz_end, child.position);
      pos += child.position;
      sum += child.position;
      count++;
      
      // if ( v == 0 ) {
      //   printf(" +++ bo nguoi ta v: %d - child: %d, c_pos: %d - %d, start: %d, end: %d\n", v, boost::target(*ei.first, graph), child.position, node.position, sz_start, sz_end);
      // }
    }

    if ( v == 0 ) {
      printf(" +++ bo nguoi ta v: %d, count: %d - %d, start: %d, end: %d\n", v, count, node.position, sz_start, sz_end);
    }
    
    pos = count == 0 ? 0 : pos/count;
    if ( pos != 0 )
      node.position = pos;
        
    if ( sz_start != sz_end && sz_start != node.position && sz_end != node.position  && count != 0 ) {
      rsOut << "chr" << chromosome << SEPARATOR
            // << v  << SEPARATOR
            // << count << SEPARATOR
            << sz_start << SEPARATOR
            << sz_end << std::endl;
    }
  }
  rsOut.close();
}

// void assureGraphPositions(Graph& graph) {  
//   for ( auto vi = boost::vertices(graph); vi.first != vi.second; ++vi.first ) {
//     auto v = *vi.first;
//     auto node = graph[v];
   
//   }
// }








