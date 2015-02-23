/****************************************************************************************
 * File: DataLoad.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 12/06/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_DATALOAD_HPP
#define SAMOGWAS_DATALOAD_HPP

#include <vector>
#include <string>
#include <string>
#include <boost/lexical_cast.hpp>
#include <iostream>

#include "utils/csv_parser.hpp"
#include "utils/matrix_utils.hpp"
namespace samogwas
{

typedef std::vector< std::vector<int> > Matrix;
typedef std::shared_ptr<Matrix> MatrixPtr;
// template< typename T >

inline MatrixPtr loadDataTable( const std::string& infile,
                         const char& sep = ',',
                         const char& quote = '"' ) {

  MatrixPtr dt(new Matrix());
  std::ifstream matrixFile(infile.c_str());
  if (!matrixFile) {
    printf("file data %s not existing\n", infile.c_str());
    exit(-1);
  }
  std::cout << "loading data" << std::endl << std::endl;

  dt->reserve(100000);

  utility::CSVIterator<int> matrixLine(matrixFile);
  
  for( ; matrixLine != utility::CSVIterator<int>(); ++matrixLine ) {         
    std::vector<int> row(matrixLine->size(), 0);
    for (unsigned i = 0; i < matrixLine->size(); ++i) {
      row[i] = matrixLine->at(i);
    }
    dt->push_back(row);    
  }

  dt->resize(dt->size());
  size_t ncols = dt->empty() ? 0 : (*dt)[0].size();
  std::cout << "done loading matrix of (" << dt->size() << "," << ncols << ")" << std::endl << std::endl;

  return dt;
}

inline void loadLabelPosition( std::vector< std::string > & labels,
                        std::vector<unsigned>& ids,
                        std::vector< int >& positions,
                        const std::string& infile )  {
  std::ifstream labPosFile(infile.c_str());
  if (!labPosFile) {
    printf("file lab-post %s not existing\n", infile.c_str());
    exit(-1);
  }
  printf("beginning loading label...\n");
  std::vector<std::string>().swap(labels); //lab2Pos.clear();
  std::vector<int>().swap(positions); //.clear();
  utility::CSVIterator<std::string> labPosLine(labPosFile);// ++labPosLine;
  int id = 0;
  for( ; labPosLine != utility::CSVIterator<std::string>(); ++labPosLine ) {
    std::string label =  (*labPosLine)[2];
    int position = boost::lexical_cast<int>( (*labPosLine)[3]);
    unsigned id = boost::lexical_cast<unsigned>( (*labPosLine)[1]);
    ids.push_back(id);
    labels.push_back(label);
    positions.push_back(position);
  }

  std::cout << "done loading. loaded " << labels.size() << " variables.\n";
}

inline void loadLabelPosition2( std::vector< std::string > & labels,
                                std::vector<unsigned>& ids,
                                std::vector< int >& positions,
                                const std::string& infile )  {
  std::ifstream labPosFile(infile.c_str());
  if (!labPosFile) {
    printf("file lab-post %s not existing\n", infile.c_str());
    exit(-1);
  }
  printf("beginning loading label now...\n");
  std::vector<std::string>().swap(labels); //lab2Pos.clear();
  std::vector<int>().swap(positions); //.clear();
  utility::CSVIterator<std::string> labPosLine(labPosFile);// ++labPosLine;
  int id = 0;
  for( ; labPosLine != utility::CSVIterator<std::string>(); ++labPosLine ) {
    // printf("id: %s, label: %s, pos: %s\n",  (*labPosLine)[0].c_str(),  (*labPosLine)[1].c_str(), 
    //        (*labPosLine)[2].c_str());
    std::string label =  (*labPosLine)[1];
    int position = boost::lexical_cast<int>( (*labPosLine)[2]);
    unsigned id = boost::lexical_cast<unsigned>( (*labPosLine)[0]);
    ids.push_back(id);
    labels.push_back(label);
    positions.push_back(position);
  }

  std::cout << "done loading. loaded " << labels.size() << " variables.\n";
}

} // namespace clusteringends here. clustering

/****************************** IMLEMENTATION BELOW THIS POINT **************************/
namespace clustering
{


} // namespace clusteringends here. clustering

/****************************************************************************************/
#endif // SAMOGWAS_DATALOAD_HPP
