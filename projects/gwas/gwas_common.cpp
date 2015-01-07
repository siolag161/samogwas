
#include <boost/program_options.hpp>

#include "gwas_common.hpp"

using namespace samogwas;

void assureGraphPositions( Graph& graph ) {
  for ( auto vi = boost::vertices(graph); vi.first != vi.second; ++vi.first ) {
    auto v = *vi.first;
    Node& node = graph[v];
    size_t sz_start = std::numeric_limits<int>::max(), sz_end = 0;
    size_t sum = 0, pos = 0, count = 0;
    for ( auto ei = boost::out_edges(v, graph); ei.first != ei.second; ++ei.first ) {
      auto child = graph[boost::target(*ei.first, graph)];
      sz_start = std::min(sz_start, (size_t)child.position);
      sz_end = std::max(sz_end, (size_t)child.position);
      pos += child.position;
      sum += child.position;
      count++;
    }
    
    pos = count == 0 ? 0 : pos/count;
    if ( pos != 0 )
      node.position = pos;
  }

  // return graph;
}



char* current_date()
{
  time_t rawtime;
  struct tm * timeinfo;
  char* buffer = new char[80];
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  strftime (buffer,80,"%Y_%m_%d_%H_%M_%S",timeinfo);

  return buffer;
}

boost::filesystem::path outputDir( std::string& outputDir, bool hasDate  ) {
  auto path = boost::filesystem::absolute(outputDir);
  if (hasDate) {
    path /= current_date();
  }
  boost::filesystem::create_directories(path);
  return path;
}


PhenoVecPtr loadPhenotype(std::string& phenoFile) {
  auto phenoVec = std::make_shared<PhenoVec>();

  std::ifstream labPosFile(phenoFile.c_str());
  if (!labPosFile) {
    std::cout << "phenotype file does not exist. program aborted\n";
    exit(-1);
  }
  utility::CSVIterator<std::string> labPosLine(labPosFile);// ++labPosLine;
  for( ; labPosLine != utility::CSVIterator<std::string>(); ++labPosLine ) {    
    int pheno = boost::lexical_cast<int>( (*labPosLine)[0]);
    phenoVec->push_back(pheno);
  }

  return phenoVec;
}


std::vector<int> getGraphParent( const samogwas::Graph& graph ) {
  std::vector<int> parent( boost::num_vertices(graph), -1 );
  for ( auto ei = boost::edges(graph); ei.first != ei.second; ++ei.first ) {
    auto source = boost::source(*ei.first, graph);
    auto target = boost::target(*ei.first, graph);
    parent[target] = source;
  }
  
  return parent;
}


////////////////
std::vector<int> countClusterSiblings( samogwas::Graph& graph ) {
  std::vector<int> sibling_count( boost::num_vertices(graph) , 0);
  
  for ( auto vi = boost::vertices(graph); vi.first != vi.second; ++vi.first ) {
    auto v = *vi.first;
    Node& node = graph[v];

    auto ei = boost::out_edges(v, graph);
    int count = std::distance(ei.first, ei.second);
    for ( ; ei.first != ei.second; ++ei.first ) {
      auto child_idx = boost::target(*ei.first, graph);
      sibling_count[child_idx] = count-1;
    }    
  }

  return sibling_count;                                  
}
