
#include "tools.hpp"

#include <string>
#include <map>
#include <fstream>
#include <boost/lexical_cast.hpp>


#include "utils/csv_parser.hpp"
#include "utils/matrix_utils.hpp"
#include "clustering/clustering.hpp"

using namespace samogwas;

Partition read_clustering_from_haploview(std::string& haploFile, std::string& mapFile) {
  std::map< std::string, int > label2id;
  std::map<int, std::string> id2label;
  std::vector<int> par2global;
  std::map<int,int> global2par;
  std::cout << "reading map...\n";
  read_map( label2id, id2label, par2global, global2par, mapFile );
  Partition partition;
  std::cout << "loading partition...\n";

  read_partition( partition, haploFile, label2id, id2label, par2global, global2par );  
  return partition;
}

void read_map( std::map<std::string, int>& label2id, 
               std::map<int, std::string>& id2label,
               std::vector<int>& par2global, std::map<int,int>& global2par, std::string mapFileName ) {
  std::ifstream mFile( mapFileName.c_str());
  if (!mFile) {
    exit(-1);
    //return;
  }

  unsigned nrows = std::count( std::istreambuf_iterator<char>(mFile), 
                               std::istreambuf_iterator<char>(), '\n' );
  mFile.seekg (0, std::ios::beg);
  
  utility::CSVIterator<std::string> mLine(mFile);// ++labPosLine;
  int par_id = 0;
  par2global.resize(nrows, 0);
  for( ; mLine != utility::CSVIterator<std::string>(); ++mLine ) {
    int global_id =  boost::lexical_cast<int>((*mLine)[1]);
    std::string label = (*mLine)[2];
    label2id[label]= global_id;
    id2label[global_id] = label;
    par2global[par_id] = global_id;
    global2par[global_id] = par_id;
    // if (label=="rs2098685") {
    //   printf("%s - local_id: %d, global_id: %d\n", label.c_str(), par_id, global_id);
    //   printf("BOOOO: %d-%s clusted\n", global_id, id2label[global_id].c_str());
    // }
        ++par_id;

  }
  mFile.close();
  std::cout << "loaded..." << label2id.size() << " vars\n";
}

void read_partition( samogwas::Partition& partition, std::string haploFile,
                     std::map< std::string, int >& label2id,
                     std::map<int, std::string>& id2label,
                     std::vector<int>& par2global, std::map<int,int>& global2par ) {
  printf("haplo file: %s\n", haploFile.c_str());
  std::ifstream hFile( haploFile.c_str());
  if (!hFile) return;

  std::vector<int> item_flag( par2global.size(), -1);
  
  utility::CSVIterator<std::string> hLine(hFile);// ++labPosLine;
  int cluster = 0;
  int clustered = 0;
  for( ; hLine != utility::CSVIterator<std::string>(); ++hLine ) {
    for ( int i = 0; i < hLine->size(); ++i ) {
      std::string label = (*hLine)[i];
      int global_id = label2id[label];
      int par_id = global2par[global_id]; 
      partition.setLabel( par_id, cluster );
      item_flag[par_id] = cluster;
      ++clustered;      
    }
    cluster++;
  }

  clustered = 0;
  for ( int i = 0; i < item_flag.size(); ++i )
  {
    int global_id = par2global[i];
    if (item_flag[i] == -1) {
      partition.setLabel( i, cluster++ );
      clustered++;
    }
  }
  hFile.close();
}

void write_partition( Partition& partition, std::map<int, std::string>& id2label,
                      std::vector<int>& par2global,
                      std::map<int,int>& global2par, std::string& outPath ) {
  std::ofstream pFile( outPath.c_str());
  int singleton = partition.nbrClusters();
  for ( size_t par_id = 0; par_id < partition.nbrClusters(); ++par_id ) {
    int clust =  partition.getLabel( par_id );        
    pFile << "\"chr02\"" << ","
          << par2global[ par_id ] << ","
          << id2label[ par2global[ par_id ] ] << ","
          << clust << ",\"haplo-plink\"" << std::endl;  
  }
  pFile.close();
}
