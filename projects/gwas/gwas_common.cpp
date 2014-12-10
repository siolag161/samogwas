
#include <boost/program_options.hpp>

#include "gwas_common.hpp"

using namespace samogwas;

void assureGraphPositions( Graph& graph ) {
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

boost::filesystem::path outputDir( std::string& outputDir ) {
  auto path = boost::filesystem::absolute(outputDir);
  path /= current_date();
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
