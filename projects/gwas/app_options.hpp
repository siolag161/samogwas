/****************************************************************************************
 * File: options.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 09/07/2014

 ***************************************************************************************/
#ifndef FLTM_OPTIONS_HPP
#define FLTM_OPTIONS_HPP

#include <string>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp> // to obtain the program's name
#include "fltm/core_fltm.hpp"
#include "utils/option_printer.hpp"
#include "utils/custom_option_desc.hpp"

namespace po = boost::program_options;   

namespace samogwas
{

struct Options {
  /* Input */
  std::string inputDataFile;
  std::string inputImputedDataFile;
  std::string inputLabelFile;
  std::string inputPheno;

  std::string graphFile;
  std::string bayesVertices;  
  std::string bayesDist;

  std::string scoreFile;
  std::string thresFile;

  double overall_thres;

  std::string outputDir;
  
};


inline Options getProgramOptions(int argc, char** argv) {
  Options result;
  std::string appName = boost::filesystem::basename(argv[0]);
  po::options_description optDesc("Options");

  try  {
    /** Define and parse the program options 
     */
    optDesc.add_options()
        ("help,h", "Print help messages")        
        ("in_dat,i", po::value<std::string>(&result.inputDataFile)->required(), "Input Data File")
        // ("in_imputed,i", po::value<std::string>(&result.inputDataFile)->required(), "Input Imputed Data File")
        ("in_pheno,p", po::value<std::string>(&result.inputPheno)->required(), "Input Pheno File")
        ("in_graph,g", po::value<std::string>(&result.graphFile)->required(), "Input Graph File")
        ("in_bayes_vertex,v", po::value<std::string>(&result.bayesVertices)->required(), "Input Bayes File")
        ("in_bayes_dist,d", po::value<std::string>(&result.bayesDist)->required(), "Input Dist File")

        ("in_lab,l", po::value<std::string>(&result.inputLabelFile)->required(), "Input Label File")
        ("score_file,s", po::value<std::string>(&result.scoreFile)->required(), "Input Score File")
        ("thres_file,t", po::value<std::string>(&result.thresFile)->required(), "Input Threshold File")
        
        ("overral_thres,r", po::value<double>(&result.overall_thres)->required(), "Input Overral Threshold")

        ("outDir,o", po::value<std::string>(&result.outputDir)->required(), "Output Dir")

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


} // namespace fltmends here. fltm

/****************************************************************************************/
#endif // FLTM_OPTIONS_HPP
