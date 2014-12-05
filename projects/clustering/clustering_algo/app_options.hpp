/****************************************************************************************
 * File: app_options.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 05/12/2014

 ***************************************************************************************/
#ifndef _APP_OPTIONS_HPP
#define _APP_OPTIONS_HPP

struct ApplicationOptions
{
  std::string dataInFile; // input filename
  std::string labPosInFile; // input filename
  std::string outputDir; // input filename

  unsigned maxDist;  
  double simi;

  std::string configFile;
  
  ApplicationOptions() {}
};

/****************************************************************************************/
#endif // _APP_OPTIONS_HPP
