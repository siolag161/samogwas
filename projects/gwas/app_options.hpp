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
  int chromosome;

  std::string graphFile;
  std::string bayesVertices;  
  std::string bayesDist;

  std::string mappingFile;
  std::string outputDir;

  int task;
  
};

struct Tool_Options: public Options {  
  std::string scoreFile;
  std::string thresFile;

  double overall_thres;
};

struct GWAS_Options: public Options {
  double threshold;
  int permutations;
};

//////////////////////////////////////////////////////////////////////////////
inline GWAS_Options getGwasProgramOptions(int argc, char** argv) {
  GWAS_Options result;
  std::string appName = boost::filesystem::basename(argv[0]);
  po::options_description optDesc("Options");

  try  {
    /** Define and parse the program options 
     */
    optDesc.add_options()
        ("help,h", "Print help messages")
        ("chr,c", po::value<int>(&result.chromosome)->default_value(2), "chromosome")

        ("in_dat,i", po::value<std::string>(&result.inputDataFile)->required(), "Input Data File")
        ("in_lab,l", po::value<std::string>(&result.inputLabelFile)->required(), "Input Label File")        

        // ("in_imputed,i", po::value<std::string>(&result.inputDataFile)->required(), "Input Imputed Data File")

        ("in_pheno,p", po::value<std::string>(&result.inputPheno)->required(), "Input Pheno File")
        ("in_graph,g", po::value<std::string>(&result.graphFile)->required(), "Input Graph File")
        ("in_bayes_vertex,v", po::value<std::string>(&result.bayesVertices)->required(), "Input Bayes File")
        ("in_bayes_dist,d", po::value<std::string>(&result.bayesDist)->required(), "Input Dist File")


        ("permutations,n", po::value<int>(&result.permutations)->default_value(1000), "Nbr Permutations")
        ("threshold,t", po::value<double>(&result.threshold)->default_value(.005), "threshold")
        ("mappingFile,m", po::value<std::string>(&result.mappingFile)->required(), "SNP - RS Mapping File")

        ("task,k", po::value<int>(&result.task)->default_value(0), "task. 0: gwas_good_parent,  1: gwas_test")
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



////////////////////////////////////////////////////////////
inline Tool_Options getToolProgramOptions(int argc, char** argv) {
  Tool_Options result;
  std::string appName = boost::filesystem::basename(argv[0]);
  po::options_description optDesc("Options");

  try  {
    /** Define and parse the program options 
     */
    optDesc.add_options()
        ("help,h", "Print help messages")        
        ("in_dat,i", po::value<std::string>(&result.inputDataFile)->required(), "Input Data File")
        // ("in_imputed,i", po::value<std::string>(&result.inputDataFile)->required(), "Input Imputed Data File")
        ("chr,c", po::value<int>(&result.chromosome)->default_value(2), "chromosome")

        ("in_pheno,p", po::value<std::string>(&result.inputPheno)->required(), "Input Pheno File")
        ("in_graph,g", po::value<std::string>(&result.graphFile)->required(), "Input Graph File")
        ("in_bayes_vertex,v", po::value<std::string>(&result.bayesVertices)->required(), "Input Bayes File")
        ("in_bayes_dist,d", po::value<std::string>(&result.bayesDist)->required(), "Input Dist File")

        ("in_lab,l", po::value<std::string>(&result.inputLabelFile)->required(), "Input Label File")
        ("score_file,s", po::value<std::string>(&result.scoreFile)->required(), "Input Score File")
        ("thres_file,t", po::value<std::string>(&result.thresFile)->required(), "Input Threshold File")
        
        ("overral_thres,r", po::value<double>(&result.overall_thres)->required(), "Input Overral Threshold")

        ("outDir,o", po::value<std::string>(&result.outputDir)->required(), "Output Dir")
        ("task,k", po::value<int>(&result.task)->default_value(0), "task. 0: gwas; 1: obtain: regions")

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
