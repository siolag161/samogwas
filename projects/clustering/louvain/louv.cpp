
#include "utils/timer_utils.hpp"
#include "utils/option_printer.hpp" //for printing
#include "utils/csv_parser.hpp"

#include <string>
#include <iostream> 
#include <fstream>
#include <cstdio>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp> // to obtain the program's nam

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

namespace utl = utility;
namespace po = boost::program_options;   
using namespace samogwas;     
using namespace samogwas::louvain;

typedef std::vector< std::vector<int> > Matrix; // We consider here only vector of int is relevant
typedef std::string Label;
typedef int Position;
typedef MutInfoSimilarity<Matrix> MutInfoSimi;

ApplicationOptions getProgramOptions(int argc, char** argv);
void saveClustering( const Partition& partition, const std::vector<unsigned>& ids, std::string clustFN );

int main(int argc, char** argv) {
  auto progOpt = getProgramOptions(argc, argv);

  utl::Timer timer, totalTimer; timer.start(); totalTimer.start();
  std::vector<Label> labels; std::vector<Position> positions; std::vector<unsigned> ids;
  std::cout << "loading data from " <<  progOpt.dataInFile << std::endl; // todo: logging
  auto matrix = loadDataTable ( progOpt.dataInFile );

  loadLabelPosition( labels, ids, positions, progOpt.labPosInFile );
  std::cout << "data loaded. rows: " << utl::nrows(*matrix) << ", columns: "
            << utl::ncols(*matrix) << ". takes: " <<  timer.display() << std::endl << std::endl; // todo: logging
  timer.restart();

  printf("Parameters: maxDist: %u, simi: %.2f\n",  progOpt.maxDist, progOpt.simi );
  std::string type_Data = progOpt.simi > 0 ? "binary" : "real";
  std::cout << "performing clustering with (" << type_Data << ")." << std::endl;

  // MutInfoSimi s(matrix, positions, progOpt.maxDist, progOpt.simi);
  std::shared_ptr<SimilarityMatrix> simi( new MutInfoSimi(matrix, positions, progOpt.maxDist, progOpt.simi) );
  // // std::shared_ptr<Graph> g(new Graph(simi));
  auto louv = std::make_shared<MethodLouvain>(simi);
  auto clustering = louv->run();
  
  boost::filesystem::path outputPath = boost::filesystem::absolute(progOpt.outputDir);
  std::string outputFileName = outputPath.string();

  boost::filesystem::create_directories(outputPath);
  char clustering_fn[256];
  char optBuf[80];
  sprintf( optBuf, "%d_%s",  progOpt.maxDist, type_Data.c_str() );
  sprintf( clustering_fn, "louvain_clustering_%s.txt", optBuf );
  outputFileName = (outputPath / clustering_fn).string();

  std::cout << "clustering done, takes " << timer.display() << std::endl;
  
  std::cout << "writing result now...\n";
  saveClustering( clustering, ids, outputFileName);
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

        // ("cast,c", po::value<double>(&result.CAST)->required(), "CAST")
        ("simi,s", po::value<double>(&result.simi)->required(), "Similarity")
        ("maxDist,m", po::value<unsigned>(&result.maxDist)->required(), "Max Distance")
        // ("verbose,v", po::value<int>(&result.verbose)->default_value(0), "Verbose")
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
