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
#include "clustering/louvain/louv.hpp"

#include "distance/dissimilarity.hpp"
#include "distance/similarity.hpp"
#include "distance/comparable.hpp"

#include "utils/logs_utils.hpp"
#include "utils/matrix_utils.hpp"

#include "app_options.hpp"
#include "data_load.hpp"
#include "clustering_parameters.hpp"

namespace utl = utility;
namespace po = boost::program_options;   
using namespace samogwas;     
using namespace samogwas::louvain;

typedef std::vector< std::vector<int> > Matrix; // We consider here only vector of int is relevant
typedef std::string Label;
typedef int Position;
typedef MutInfoSimilarity<Matrix> MutInfoSimi;

ApplicationOptions getProgramOptions(int argc, char** argv);

std::string clusteringFileName( ClustAlgoPtr algo, ApplicationOptions& progOpt, boost::filesystem::path& path );
std::string statFileName( ApplicationOptions& progOpt, boost::filesystem::path& path );
boost::filesystem::path outputDir( ApplicationOptions& progOpt );
char* current_date();
void saveClustering( const Partition& partition, const std::vector<unsigned>& ids, std::string clustFN );

void singletonClustering( const Partition& partition );

int main(int argc, char** argv) {
  auto progOpt = getProgramOptions(argc, argv);

  utl::Timer timer, totalTimer; timer.start(); totalTimer.start();
  std::vector<Label> labels; std::vector<Position> positions; std::vector<unsigned> ids;
  std::cout << "loading data from " <<  progOpt.dataInFile << std::endl; // todo: logging
  auto matrix = loadDataTable ( progOpt.dataInFile );

  loadLabelPosition( labels, ids, positions, progOpt.labPosInFile );
  std::cout << "data loaded. rows: " << utl::nrows(*matrix) << ", columns: "
            << utl::ncols(*matrix) << ". takes: " <<  timer.display() << std::endl << std::endl; // todo: logging

  printf("Parameters: maxDist: %u, simi: %.2f\n",  progOpt.maxDist, progOpt.simi );
  std::string type_Data = progOpt.simi > 0 ? "binary" : "real";
  std::cout << "performing clustering with (" << type_Data << ")." << std::endl;

  std::ifstream cfgIs(progOpt.configFile);
  printf("reading algo...\n");
  auto clust_algos = read_clustering_algos(matrix, positions, progOpt.maxDist, progOpt.simi, cfgIs);
  printf("done reading algo...\n");

  cfgIs.close();

  auto outPath = outputDir(progOpt);
  std::ofstream statIs(statFileName(progOpt,outPath));


  statIs << "clustering: " << labels.size() << " variables\n";
  for ( auto algo: clust_algos ) {
    std::cout << "\n-----------------------------------------------\n\n";

    timer.restart();
    printf("starting to perform %s now ", algo->name());
    statIs << "clustering: " << algo->name();
    
    auto clustering = algo->run();
    printf("done %s. got %d in %s\n ", algo->name(), clustering.nbrClusters(), timer.display().c_str()) ;

    // std::cout << "clustering %s done, takes " << timer.display() << std::endl;
    statIs << "clustering done. produces: " << clustering.nbrClusters() << "  clusters and takes " << timer.display() << std::endl;
    singletonClustering(clustering);
    std::cout << std::endl;
    auto clt_fn = clusteringFileName( algo, progOpt, outPath );
    std::cout << "writing result now...\n";
    saveClustering( clustering, ids, clt_fn );
    std::cout << "\n-----------------------------------------------\n\n";
  }

  statIs.close();
  cfgIs.close();

  // std::shared_ptr<SimilarityMatrix> simi( new MutInfoSimi(matrix, positions, progOpt.maxDist, progOpt.simi) );

  // auto louv = std::make_shared<MethodLouvain>(simi);
  // auto clustering = louv->run();
  
  // boost::filesystem::path outputPath = boost::filesystem::absolute(progOpt.outputDir);
  // std::string outputFileName = outputPath.string();

  // boost::filesystem::create_directories(outputPath);
  // char clustering_fn[256];
  // char optBuf[80];
  // sprintf( optBuf, "%d_%s",  progOpt.maxDist, type_Data.c_str() );
  // sprintf( clustering_fn, "louvain_clustering_%s.txt", optBuf );
  // outputFileName = (outputPath / clustering_fn).string();

  // std::cout << "clustering done, takes " << timer.display() << std::endl;
  
  // std::cout << "writing result now...\n";
  // saveClustering( clustering, ids, outputFileName);
  // std::cout << "takes: " << timer.display() << std::endl;
}

/////////////////
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
        ("outputDir,o", po::value<std::string>(&result.outputDir)->required(), "Output Directory")

        ("simi,s", po::value<double>(&result.simi)->required(), "Similarity")
        ("maxDist,m", po::value<unsigned>(&result.maxDist)->required(), "Max Distance")
        ("configFile,c", po::value<std::string>(&result.configFile)->required(), "configFile")
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

static const std::string ID = "id"; static const std::string LABEL = "label";
static const std::string LATENT = "latent"; static const std::string PARENT = "parent";
static const std::string LEVEL = "level"; static const std::string POSITION = "position";
static const std::string CARDINALITY = "cardinality";
static const std::string PARENT_ID = "parent_id";
static const char SEPARATOR = ',';


void saveClustering( const Partition& partition, const std::vector<unsigned>& ids, std::string clustFN ) {  
  std::ofstream clustOut(clustFN);
  clustOut << ID << SEPARATOR << LATENT << SEPARATOR << PARENT << SEPARATOR
           << LEVEL << SEPARATOR << CARDINALITY << "\n";  // writes header

  std::cout << "saving clustering of " << partition.nbrClusters() << " clusters into " << clustFN << std::endl;
  unsigned max_id = ids[ ids.size() - 1 ] + 1;

  for ( size_t var = 0; var < partition.nbrItems(); ++var ) {
    // if ( dbscan.get_labels()[var] >= 0)
    auto lab = partition.getLabel(var);
    clustOut << ids[var] << SEPARATOR << 0 << SEPARATOR << ( lab + max_id ) << SEPARATOR
             << 0 << SEPARATOR << 3 << "\n";
    // else
    //   clustOut << ids[var] << SEPARATOR << 0 << SEPARATOR << ( dbscan.get_labels()[var] ) << SEPARATOR
    //            << 0 << SEPARATOR << 3 << "\n";
  }

  for ( size_t var = max_id; var < max_id + partition.nbrClusters(); ++var ) {
    clustOut << var << SEPARATOR << 1 << SEPARATOR << -1 << SEPARATOR
             << 1 << SEPARATOR << 3 << "\n";    
  }

  clustOut.close();

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

boost::filesystem::path outputDir( ApplicationOptions& progOpt ) {
  auto path = boost::filesystem::absolute(progOpt.outputDir);
  // boost::filesystem::create_directories(path);
  path /= current_date();
  boost::filesystem::create_directories(path);

  return path;
}


std::string statFileName( ApplicationOptions& progOpt, boost::filesystem::path& path ) {
  std::string statFileName = path.string();
  std::string type_data = progOpt.simi > 0 ? "binary" : "real";
    
  char clustering_fn[256];
  char optBuf[80];
  sprintf( optBuf, "%d_%s",  progOpt.maxDist, type_data.c_str() );
  statFileName = (path / optBuf).string();

  return statFileName;
}


std::string clusteringFileName( ClustAlgoPtr algo, ApplicationOptions& progOpt, boost::filesystem::path& path ) {
  std::string outputFileName = path.string();
  std::string type_data = progOpt.simi > 0 ? "binary" : "real";
    
  char clustering_fn[256];
  char optBuf[80];
  sprintf( optBuf, "%d_%s",  progOpt.maxDist, type_data.c_str() );
  sprintf( clustering_fn, "%s_%s.txt", algo->name(), optBuf );
  outputFileName = (path / clustering_fn).string();

  return outputFileName;
}


void singletonClustering( const Partition& partition ) {
  auto clustering = partition.to_clustering();

  int nbr_non_singletons = 0;
  int total = 0;
  double ave_count = 0.0;
  int total_elems = 0;
  for ( auto clt: clustering ) {
    total++;
    total_elems += clt.size();
    if (clt.size() > 1)
    {
      nbr_non_singletons++;
      ave_count += clt.size();
    }
  }
  ave_count /= nbr_non_singletons;
  printf("nbr non-singeleton: %d (over: %d - %d), average: %f over %d\n",
         nbr_non_singletons, total, partition.nbrClusters(), ave_count, total_elems);
}
