/****************************************************************************************
 * File: main.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 09/01/2015

 ***************************************************************************************/
#ifndef PDT_MAIN_HPP
#define PDT_MAIN_HPP


struct ApplicationOptions
{
  std::string dataInFile; // input filename
  std::string labPosInFile; // input filename

  std::string refClusteringPath;
  
  std::string outputDir; // input filename

  unsigned maxDist;  
  // double simi;

  std::string configFile;
  
  ApplicationOptions() {}
};


ApplicationOptions getProgramOptions(int argc, char** argv);


/****************************************************************************************/
#endif // PDT_MAIN_HPP
