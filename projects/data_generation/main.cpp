
#include "utils/timer_utils.hpp"
#include "utils/option_printer.hpp" //for printing
#include "utils/csv_parser.hpp"

#include <string>
#include <iostream> 
#include <fstream>
#include <cstdio>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp> // to obtain the program's nam

#include <boost/lexical_cast.hpp>

#include <random>

#include "utils/logs_utils.hpp"
#include "utils/matrix_utils.hpp"

#include "app_options.hpp"
#include "data_generation.hpp"

namespace utl = utility;
namespace po = boost::program_options;   
using namespace samogwas;     

typedef std::vector< std::vector<int> > Matrix; // We consider here only vector of int is relevant
typedef std::string Label;
typedef int Position;

ApplicationOptions getProgramOptions(int argc, char** argv);
void saveData( const Matrix& mat, const std::string output_path );
void saveLabel( const size_t nbr_elems, const std::string output_path);

/**
 *
 */
int main(int argc, char** argv) {
  auto progOpt = getProgramOptions(argc, argv);
  const int CARD = 3;
  // const int NCOLS = 5000;  
  auto data = data_gen::GenerateClusteredData( progOpt.nbrClust, progOpt.clustSize, CARD, progOpt.nbrAttrs )();

  auto outputPath = boost::filesystem::absolute(progOpt.outputDir);
  boost::filesystem::create_directories(outputPath);
  char data_fn[256], label_fn[256];
  char optBuf[80];

  sprintf( optBuf, "%d_%d_%d",  progOpt.maxDist, progOpt.nbrClust, progOpt.clustSize );
  sprintf( data_fn, "louvain_clustering_%s_data.txt", optBuf );
  sprintf( label_fn, "louvain_clustering_%s_label.txt", optBuf );

  auto data_fn_str = (outputPath / data_fn).string();
  auto label_fn_str = (outputPath / label_fn).string();
  saveLabel( data->size(), label_fn_str );
  saveData( *data, data_fn_str );
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
        // ("dinput,i", po::value<std::string>(&result.dataInFile)->required(), "Data Input filename")
        // ("lpinput,l", po::value<std::string>(&result.labPosInFile)->required(), "Label-Pos Input filename")
        ("outputDir,o", po::value<std::string>(&result.outputDir)->required(), "Output Directory")

        ("nbrAttrs,c", po::value<unsigned>(&result.nbrAttrs)->required(), "Number of attributes")

        ("maxDist,m", po::value<unsigned>(&result.maxDist)->required(), "Max distance allowed")
        ("nbrClust,n", po::value<unsigned>(&result.nbrClust)->required(), "Number of clusters")
        ("clustSize,s", po::value<unsigned>(&result.clustSize)->required(), "Number of item per cluster")

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

void saveData( const Matrix& mat, const std::string dataPath ) {
  std::cout << "saving data...\n" << std::endl;  
  ////////////////////////////////
  std::ofstream matOut(dataPath);

  for(size_t row = 0; row < utility::nrows(mat); row++) {
    for(size_t col = 0; col < utility::ncols(mat) - 1; col++) {
      matOut << mat[row][col] << ",";
    }
    if ( utility::ncols(mat) > 0) {
      matOut << mat[row][utility::ncols(mat)-1];
    }
    matOut << std::endl;
  }
  matOut.close();
}

void saveLabel( const size_t nbr_elems, const std::string output_path) {
  std::ofstream labPosOut(output_path);
  // vertex_iterator vi, vi_end;
  // int latId = input.indexes[input.indexes.size()-1];
  // for ( boost::tie(vi, vi_end) = boost::vertices(input.graph); vi != vi_end; ++vi ) {
  //   vertex_t vertex = *vi;
  //   if ( vertex < input.indexes.size() )
  //     labPosOut << input.indexes[vertex] << "," << input.graph[vertex].label << "," << input.graph[vertex].position
  //               << "," << input.graph[vertex].variable.cardinality() << std::endl;
  //   else {
  //     labPosOut << ++latId << "," << "\"imputed-" + input.graph[vertex].label << "\"," << input.graph[vertex].position
  //               << "," << input.graph[vertex].variable.cardinality() << std::endl;
  //   }
  // }

  std::random_device rd;
  std::mt19937 gen(rd());
  int LOWER = 10000, UPPER = 15000;
  std::uniform_int_distribution<> dist(LOWER, UPPER);

  int curr =  dist(gen);
  
  for ( int i = 0; i < nbr_elems; ++ i ) {
    char label[80];
    sprintf( label, "\"lab-%d\"", i);
    std::string lab = boost::lexical_cast<std::string>(i);
    labPosOut << "2" << "," << i << "," << label << "," << curr << std::endl;
    curr += dist(gen);
  }

  labPosOut.close();
}
