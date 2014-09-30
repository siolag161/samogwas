
// #define BOOST_TEST_DYN_LINK
// #ifdef STAND_ALONE
// #   define BOOST_TEST_MODULE
// #endif  
// #include <boost/test/unit_test.hpp>
// #include <fstream>
// #include <map>
// #include <boost/random.hpp>
// #include <boost/random/normal_distribution.hpp>
// #include <math.h>
// #include <boost/lexical_cast.hpp>

// #include "CAST.hpp"
// #include "DBSCAN.hpp"
// #include "Simi.hpp"
// #include "EM.hpp"
// #include "DataMatrix.hpp"
// #include "Graph.hpp"
// #include "GraphIO.hpp"

// using namespace fltm;
// using namespace boost;
// class Data   
// { 
// };

// BOOST_FIXTURE_TEST_SUITE( Test_GraphLoad, Data ) 

// BOOST_AUTO_TEST_CASE( Test_Gaussian_K_DBSCAN_10 ) {
//   std::cout << "TESTING................\n\n";
//   BayesGraphLoad loader;
//   Graph graph;
//   loader( graph, "/home/pdt/BIOSTEC/rs/nv/fltm_CAST_0.500_imputedLab.out",
//           "/home/pdt/BIOSTEC/rs/nv/fltm_CAST_0.500_bayes_vertex.out",
//           "/home/pdt/BIOSTEC/rs/nv/fltm_CAST_0.500_bayes_dist.out");

//   typedef std::map<std::string, int> LabPosMap;
//   LabPosMap lpMap;
//   LabPosMap id2Pos;

//   std::ifstream labPosFile("/home/pdt/BIOSTEC/rs/nv/fltm_CAST_0.500_imputedLab.out");    
//   utility::CSVIterator<std::string> labPosLine(labPosFile); // ++labPosLine;
//   for( ; labPosLine != utility::CSVIterator<std::string>(); ++labPosLine ) {
//     size_t id = boost::lexical_cast<size_t>( (*labPosLine)[LP_ID] );
//     std::string label = (*labPosLine)[LP_LABEL] ;
//     int position = boost::lexical_cast<int>( (*labPosLine)[LP_POSITION] );
//     lpMap[label] = id;
//     id2Pos[label] = position;
//   } 
//   std::ofstream pos("/home/pdt/Desktop/CAST_real_pos.csv");
//   for ( int i = 0; i < boost::num_vertices(graph); ++i) {
//     Node& n = graph[i];
//     if ( n.level == 1 ) {
//       auto col = n.jointDistribution.get_computable_object_list();    
//       int nbrChildren = col.nbrVariables() - 1;
//       if (nbrChildren > 0) {    
//         unsigned totalPos = 0;
//         for (unsigned child = 1; child < col.nbrVariables(); ++ child) {
//           LabelType label = col[child].get_variables()[0].name();
//           IndexType index = lpMap[label];
//           totalPos += graph[index].position;

//         }    
//         int position = totalPos / nbrChildren;
//         n.position = position;
//       }      
//     }
//   }

//   for ( int i = 0; i < boost::num_vertices(graph); ++i) {
//     Node& n = graph[i];
//     if ( n.level == 2 ) {
//       auto col = n.jointDistribution.get_computable_object_list();    
//       int nbrChildren = col.nbrVariables() - 1;
//       if (nbrChildren > 0) {    
//         unsigned totalPos = 0;
//         for (unsigned child = 1; child < col.nbrVariables(); ++ child) {
//           LabelType label = col[child].get_variables()[0].name();
//           IndexType index = lpMap[label];
//           totalPos += graph[index].position;
//         }    
//         int position = totalPos / nbrChildren;
//         n.position = position;
//       }      
//     }
//   }


//   for ( int i = 0; i < boost::num_vertices(graph); ++i) {
//     Node& n = graph[i];
//     if ( !n.isLeaf ) {
//       auto col = n.jointDistribution.get_computable_object_list();    
//       int nbrChildren = col.nbrVariables() - 1;
//       int maxChildrenLevel = 0;

//       if (nbrChildren > 0) {    
//         unsigned totalPos = 0;
//         for (unsigned child = 1; child < col.nbrVariables(); ++ child) {
//           LabelType label = col[child].get_variables()[0].name();
//           IndexType index = lpMap[label];
//           totalPos += graph[index].position;
//           if (maxChildrenLevel < graph[child].level) {
//             maxChildrenLevel = graph[child].level;
//           }
//         }    
//         int position = totalPos / nbrChildren;
//         n.position = position;
//         n.level = maxChildrenLevel + 1;
//       }      
//     }
//     pos << i << "," << n.position << "," <<  n.level << std::endl;
//   }
//   pos.close();  
//   std::cout << "TESTING................\n\n";
// }


// BOOST_AUTO_TEST_CASE( Test_Gaussian_K_DBSCAN_101 ) {
//   std::cout << "TESTING................\n\n";
//   BayesGraphLoad loader;
//   Graph graph;
//   loader( graph, "/home/pdt/BIOSTEC/rs/nv/fltm_DBSCAN_2_0.300_imputedLab.out",
//           "/home/pdt/BIOSTEC/rs/nv/fltm_DBSCAN_2_0.300_bayes_vertex.out",
//           "/home/pdt/BIOSTEC/rs/nv/fltm_DBSCAN_2_0.300_bayes_dist.out");

//   typedef std::map<std::string, int> LabPosMap;
//   LabPosMap lpMap;
//   LabPosMap id2Pos;

//   std::ifstream labPosFile("/home/pdt/BIOSTEC/rs/nv/fltm_DBSCAN_2_0.300_imputedLab.out");    
//   utility::CSVIterator<std::string> labPosLine(labPosFile); // ++labPosLine;
//   for( ; labPosLine != utility::CSVIterator<std::string>(); ++labPosLine ) {
//     size_t id = boost::lexical_cast<size_t>( (*labPosLine)[LP_ID] );
//     std::string label = (*labPosLine)[LP_LABEL] ;
//     int position = boost::lexical_cast<int>( (*labPosLine)[LP_POSITION] );
//     lpMap[label] = id;
//     id2Pos[label] = position;
//   } 
//   std::ofstream pos("/home/pdt/Desktop/DBSCAN_pos.csv");
//   for ( int i = 0; i < boost::num_vertices(graph); ++i) {
//     Node& n = graph[i];
//     if ( n.level == 1 ) {
//       auto col = n.jointDistribution.get_computable_object_list();    
//       int nbrChildren = col.nbrVariables() - 1;
//       if (nbrChildren > 0) {    
//         unsigned totalPos = 0;
//         for (unsigned child = 1; child < col.nbrVariables(); ++ child) {
//           LabelType label = col[child].get_variables()[0].name();
//           IndexType index = lpMap[label];
//           totalPos += graph[index].position;

//         }    
//         int position = totalPos / nbrChildren;
//         n.position = position;
//       }      
//     }
//   }

//   for ( int i = 0; i < boost::num_vertices(graph); ++i) {
//     Node& n = graph[i];
//     if ( n.level == 2 ) {
//       auto col = n.jointDistribution.get_computable_object_list();    
//       int nbrChildren = col.nbrVariables() - 1;
//       if (nbrChildren > 0) {    
//         unsigned totalPos = 0;
//         for (unsigned child = 1; child < col.nbrVariables(); ++ child) {
//           LabelType label = col[child].get_variables()[0].name();
//           IndexType index = lpMap[label];
//           totalPos += graph[index].position;
//         }    
//         int position = totalPos / nbrChildren;
//         n.position = position;
//       }      
//     }
//   }


//   for ( int i = 0; i < boost::num_vertices(graph); ++i) {
//     Node& n = graph[i];
//     if ( !n.isLeaf ) {
//       auto col = n.jointDistribution.get_computable_object_list();    
//       int nbrChildren = col.nbrVariables() - 1;
//       int maxChildrenLevel = 0;

//       if (nbrChildren > 0) {    
//         unsigned totalPos = 0;
//         for (unsigned child = 1; child < col.nbrVariables(); ++ child) {
//           LabelType label = col[child].get_variables()[0].name();
//           IndexType index = lpMap[label];
//           totalPos += graph[index].position;
//           if (maxChildrenLevel < graph[child].level) {
//             maxChildrenLevel = graph[child].level;
//           }
//         }    
//         int position = totalPos / nbrChildren;
//         n.position = position;
//         n.level = maxChildrenLevel + 1;
//       }      
//     }
//     pos << i << "," << n.position << "," <<  n.level << std::endl;
//   }
//   pos.close();  
//   std::cout << "TESTING................\n\n";
// }


// BOOST_AUTO_TEST_CASE( Test_Gaussian_K_DBSCAN_102 ) {
//   std::cout << "TESTING................\n\n";
//   BayesGraphLoad loader;
//   Graph graph;
//   loader( graph, "/home/pdt/Desktop/fltm/CAST_bin/fltm_CAST_0.500_imputedLab.out",
//           "/home/pdt/Desktop/fltm/CAST_bin/fltm_CAST_0.500_bayes_vertex.out",
//           "/home/pdt/Desktop/fltm/CAST_bin/fltm_CAST_0.500_bayes_dist.out");

//   typedef std::map<std::string, int> LabPosMap;
//   LabPosMap lpMap;
//   LabPosMap id2Pos;

//   std::ifstream labPosFile("/home/pdt/Desktop/fltm/CAST_bin/fltm_CAST_0.500_imputedLab.out");    
//   utility::CSVIterator<std::string> labPosLine(labPosFile); // ++labPosLine;
//   for( ; labPosLine != utility::CSVIterator<std::string>(); ++labPosLine ) {
//     size_t id = boost::lexical_cast<size_t>( (*labPosLine)[LP_ID] );
//     std::string label = (*labPosLine)[LP_LABEL] ;
//     int position = boost::lexical_cast<int>( (*labPosLine)[LP_POSITION] );
//     lpMap[label] = id;
//     id2Pos[label] = position;
//   } 
//   std::ofstream pos("/home/pdt/Desktop/CAST_bin_pos.csv");
//   for ( int i = 0; i < boost::num_vertices(graph); ++i) {
//     Node& n = graph[i];
//     if ( n.level == 1 ) {
//       auto col = n.jointDistribution.get_computable_object_list();    
//       int nbrChildren = col.nbrVariables() - 1;
//       if (nbrChildren > 0) {    
//         unsigned totalPos = 0;
//         for (unsigned child = 1; child < col.nbrVariables(); ++ child) {
//           LabelType label = col[child].get_variables()[0].name();
//           IndexType index = lpMap[label];
//           totalPos += graph[index].position;

//         }    
//         int position = totalPos / nbrChildren;
//         n.position = position;
//       }      
//     }
//   }

//   for ( int i = 0; i < boost::num_vertices(graph); ++i) {
//     Node& n = graph[i];
//     if ( n.level == 2 ) {
//       auto col = n.jointDistribution.get_computable_object_list();    
//       int nbrChildren = col.nbrVariables() - 1;
//       if (nbrChildren > 0) {    
//         unsigned totalPos = 0;
//         for (unsigned child = 1; child < col.nbrVariables(); ++ child) {
//           LabelType label = col[child].get_variables()[0].name();
//           IndexType index = lpMap[label];
//           totalPos += graph[index].position;
//         }    
//         int position = totalPos / nbrChildren;
//         n.position = position;
//       }      
//     }
//   }


//   for ( int i = 0; i < boost::num_vertices(graph); ++i) {
//     Node& n = graph[i];
//     if ( !n.isLeaf ) {
//       auto col = n.jointDistribution.get_computable_object_list();    
//       int nbrChildren = col.nbrVariables() - 1;
//       int maxChildrenLevel = 0;

//       if (nbrChildren > 0) {    
//         unsigned totalPos = 0;
//         for (unsigned child = 1; child < col.nbrVariables(); ++ child) {
//           LabelType label = col[child].get_variables()[0].name();
//           IndexType index = lpMap[label];
//           totalPos += graph[index].position;
//           if (maxChildrenLevel < graph[child].level) {
//             maxChildrenLevel = graph[child].level;
//           }
//         }    
//         int position = totalPos / nbrChildren;
//         n.position = position;
//         n.level = maxChildrenLevel + 1;
//       }      
//     }
//     pos << i << "," << n.position << "," <<  n.level << std::endl;
//   }
//   pos.close();  
//   std::cout << "TESTING................\n\n";
// }



// BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here
