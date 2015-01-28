/****************************************************************************************
 * File: tools.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 09/01/2015

 ***************************************************************************************/
#ifndef _TOOLS_HPP
#define _TOOLS_HPP

#include <string>
#include <map>
#include <vector>
#include "clustering/clustering.hpp"

void read_map( std::map<std::string, int>& label2id, 
               std::map<int, std::string>& id2label,
               std::vector<int>& par2global, std::map<int,int>& global2par, std::string mapFileName );

void read_partition( samogwas::Partition&, std::string haplo,
                     std::map< std::string, int >& label2id,
                     std::map<int, std::string>& id2label,

                     std::vector<int>& par2global, std::map<int,int>& global2par  );

void write_partition( samogwas::Partition& partition, std::map<int, std::string>& id2label,
                      std::vector<int>& par2global,
                      std::map<int,int>& global2par, std::string& outPath );


samogwas::Partition read_clustering_from_haploview(std::string& infile, std::string& mapfile);

void output_stats( samogwas::Partition&, std::string& file );

/****************************************************************************************/
#endif // _TOOLS_HPP
