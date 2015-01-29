#include "utils/timer_utils.hpp"
#include "utils/option_printer.hpp" //for printing
#include "utils/csv_parser.hpp"

#include <string>
#include <iostream> 
#include <fstream>
#include <cstdio>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp> // to obtain the program's nam
#include <ctime>
#include <iomanip> 

#include "clustering/clustering.hpp"
#include "clustering/cast.hpp"
#include "clustering/compare_measure.hpp"
#include "clustering/rand.hpp"
#include "clustering/fowlkes.hpp"
#include "clustering/clustering_quality.hpp"

#include "clustering/louvain/louv.hpp"

#include "distance/dissimilarity.hpp"
#include "distance/similarity.hpp"
#include "distance/comparable.hpp"

#include "utils/logs_utils.hpp"
#include "utils/matrix_utils.hpp"

#include "clustering_parameters.hpp"
#include "tools.hpp"
#include "main.hpp"
#include "data_load.hpp"

namespace utl = utility;
namespace po = boost::program_options;   
using namespace samogwas;     
using namespace samogwas::louvain;

typedef std::vector< std::vector<int> > Matrix; // We consider here only vector of int is relevant
typedef std::string Label;
typedef int Position;
typedef MutInfoSimilarity<Matrix> MutInfoSimi;

/////////////////////////////////////////////////////////////////////////
double averageSize(const Clustering& clustering);
int nbrSingletons(const Clustering& clustering);
inline void outputStatistics(const Partition& clustering, const std::string& algo_params, const std::string& file);
std::string getOutFilePath(std::string path, std::string inf);


////////////////////////////////////////////
int main(int argc, char** argv) {
  utl::Timer timer, totalTimer; timer.start(); totalTimer.start();

  auto progOpt = getProgramOptions(argc, argv);

  std::vector<Label> labels; std::vector<Position> positions; std::vector<unsigned> ids;
  std::cout << "loading data from " <<  progOpt.dataInFile << std::endl; // todo: logging
  auto matrix = loadDataTable ( progOpt.dataInFile );

  loadLabelPosition( labels, ids, positions, progOpt.labPosInFile );
  std::cout << "data loaded. rows: " << utl::nrows(*matrix) << ", columns: "
            << utl::ncols(*matrix) << ". takes: " <<  timer.display() << std::endl << std::endl; // todo: logging

  printf("Parameters - : - maxDist: %u, simi: %f\n",  progOpt.maxDist, -1.0 );

  boost::filesystem::path inputPath(progOpt.dataInFile);

  ClustAlgoPtr algo;
  // // // // find the best cast config
  double cast = 0.05;
  auto simi = std::make_shared<MutInfoSimi>( matrix, positions, progOpt.maxDist, -1 );      

  std::string outFilePath = getOutFilePath(progOpt.outputDir, inputPath.filename().string());

    
  std::ofstream outFile(outFilePath);
  double score = -1.0;
  while (cast < 0.9) {
    // algo = std::make_shared<CAST_Algo>(-1, cast);   
    cast += 0.05;
    algo = std::make_shared<CAST_Algo>(simi, cast);
    samogwas::Partition partition = algo->run();

    auto nbrClusters = partition.nbrClusters();
    auto clustering = partition.to_clustering();
    auto nbrSings = nbrSingletons(clustering);
    outFile << algo->name() << ";" << ""
            << nbrSings  << ";"
            << nbrClusters << ";"
            << (double)nbrSings / nbrClusters << ";"
            << averageSize(clustering) << ";"
            << std::endl;
  }

  auto diss = std::make_shared<MutInfoDiss>( matrix, positions, progOpt.maxDist, -1 );      


  for ( int minPts = 1; minPts < 4; ++minPts ) {
    for ( double eps = 0.2; eps < 0.75; ++eps ) {
      algo = std::make_shared<DBSCAN_Algo>( diss, minPts, eps );      
      samogwas::Partition partition = algo->run();
      auto nbrClusters = partition.nbrClusters();
      auto clustering = partition.to_clustering();
      auto nbrSings = nbrSingletons(clustering);
      outFile << algo->name() << ";" << ""
              << nbrSings  << ";"
              << nbrClusters << ";"
              << (double)nbrSings / nbrClusters << ";"
              << averageSize(clustering) << ";"
              << std::endl;
    }
  }
  
  outFile.close();
  // find the best dbscan config
  
  
  // find the best louvain config
  

  // write the thing here
  

  // write the thing here

  
}

ApplicationOptions getProgramOptions(int argc, char** argv)
{
  ApplicationOptions result;
  std::string appName = boost::filesystem::basename(argv[0]);
  po::options_description optDesc("Options");
      
  try  {
    /** Define and parse the program options 
     */
    optDesc.add_options()
        ("help,h", "Print help messages")
        ("dinput,i", po::value<std::string>(&result.dataInFile)->required(), "Data Input filename")
        ("lpinput,l", po::value<std::string>(&result.labPosInFile)->required(), "Label-Pos Input filename")
        // ("refpath,r", po::value<std::string>(&result.refClusteringPath)->required(), "refClusteringPath")

        ("maxDist,m", po::value<unsigned>(&result.maxDist)->required(), "Max Distance")

        ("outputDir,o", po::value<std::string>(&result.outputDir)->required(), "Output Directory")
 

        ;
    po::variables_map vm; 
    try { 
      po::store(po::command_line_parser(argc, argv).options(optDesc).run(), vm); // throws on error
      if (vm.count("help") ) {
        samogwas::OptionPrinter::printStandardAppDesc(appName,std::cout, optDesc, NULL);
        exit(1);
      }
      po::notify(vm);   	    

    } 
    catch(boost::program_options::required_option& e) /** missing arguments **/
    {
      samogwas::OptionPrinter::formatRequiredOptionError(e);
      std::cout << e.what() << std::endl << std::endl;
      samogwas::OptionPrinter::printStandardAppDesc( appName,std::cout,
                                                     optDesc, NULL);
      exit(-1);
    }

  }
  catch(std::exception e)    
  {
    std::cout << "Unhandled Exception reached the top of main: "
              << e.what() << ", application will now exit" << std::endl;

    samogwas::OptionPrinter::printStandardAppDesc(appName, std::cout, optDesc, NULL);
    exit(-1);
  }

  return result;
}

double averageSize(const Clustering& clustering) {
  auto nbrClusters = clustering.size();
  size_t nbrItems = 0;
  for (auto& clt: clustering) {
    nbrItems += clt.size();
  }

  double averageSize = (double)nbrItems/nbrClusters;
  return averageSize;
}

int nbrSingletons(const Clustering& clustering) {
  int nbrSingletons = 0;
  for (auto& clt: clustering) {
    if (clt.size()==1) {
      ++nbrSingletons;
    }
  }

  return nbrSingletons;
}

// void outputStatistics(const Clustering& clustering, const std::string& algo_params, const std::string& file) {
//   auto nbrClusters = clustering.size();

//   int nbrSingleton = 0;
//   size_t nbrItems = 0;
//   for (auto& clt: clustering) {
//     if (clt.size()==1) {
//       ++nbrSingleton;
//     }
//     nbrItems += clt.size();
//   }


//   std::ofstream outFile(file);
//   file << algo_params << std::endl;
//   file << "averageItems"

//   outFile.close();  
// }


std::string getOutFilePath(std::string path, std::string inf) {
  boost::filesystem::path outputPath = boost::filesystem::absolute(path);
  std::string outputFileName = outputPath.string();

  boost::filesystem::create_directories(outputPath);
  char fn[256];
  sprintf( fn, "statistic_%s", inf.c_str() );
  outputFileName = (outputPath / fn).string();
  std::cout << "filename: " << outputFileName << std::endl;
  
  return outputFileName;

}
