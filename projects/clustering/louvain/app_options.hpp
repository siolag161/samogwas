/****************************************************************************************
 * File: app_options.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 01/12/2014

 ***************************************************************************************/
#ifndef LOUVAIN_APP_OPTIONS_HPP
#define LOUVAIN_APP_OPTIONS_HPP

struct ApplicationOptions
{
  std::string dataInFile; // input filename
  std::string labPosInFile; // input filename
  std::string outputDir; // input filename

  unsigned maxDist;
  
  // int hasHeader; // indicates whether the input file contains a header row
  // int verbose;
  double simi;
  ApplicationOptions() {}
};

/****************************************************************************************/
#endif // LOUVAIN_APP_OPTIONS_HPP
