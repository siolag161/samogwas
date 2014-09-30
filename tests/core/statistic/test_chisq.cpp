#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
 #   define BOOST_TEST_MODULE
 #endif  
 #include <boost/test/unit_test.hpp>

 #include "statistics/test_statistics.hpp"
 #include <omp.h>
 #include <chrono>
 // #include "DataLoad.hpp"

 #include <boost/filesystem.hpp> // to obtain the program's name
 class Data 
 { 
 };
 
BOOST_FIXTURE_TEST_SUITE( Test_G2, Data ) 


// // OOST_AUTO_TEST_CASE( Test_ChiSqrt )
// // {
// //   std::vector<int> geno;
// //   for ( int i = 0; i < 12; ++i ) geno.push_back(0);
// //   for ( int i = 0; i < 5; ++i ) geno.push_back(1);
// //   for ( int i = 0; i < 7; ++i ) geno.push_back(0);
// //   for ( int i = 0; i < 7; ++i ) geno.push_back(1);
// //   std::vector<int> pheno;
// //   for ( int i = 0; i < 17; ++i ) pheno.push_back(0);
// //   for ( int i = 0; i < 14; ++i ) pheno.push_back(1);
  
// //   stats::StatisticTest<stats::CHISQ> chisq;
// //   // BOOST_CHECK_CLOSE( chisq(geno, pheno, 2, 2), 0.2415293, 0.0001 );
// // }

BOOST_AUTO_TEST_CASE(  Test_ChiSqrt_Corrected )
{
  std::vector<int> geno;
  for ( int i = 0; i < 12; ++i ) geno.push_back(0);
  for ( int i = 0; i < 5; ++i ) geno.push_back(1);
  for ( int i = 0; i < 7; ++i ) geno.push_back(0);
  for ( int i = 0; i < 7; ++i ) geno.push_back(1);
  
  std::vector<int> pheno;
  for ( int i = 0; i < 17; ++i ) pheno.push_back(0);
  for ( int i = 0; i < 14; ++i ) pheno.push_back(1);
  
  stats::StatisticTest<stats::CHISQ_YATES> chisq;
  BOOST_CHECK_CLOSE( chisq(geno, pheno, 2, 2), 0.4233054, 0.0001 );
}

// /** checks if input exists and exists on giving the error message
//  *
//  */
// void checkInputFiles( std::string path, std::string filename ) {
//   if ( !boost::filesystem::exists( path ) )
//   {
//     std::cout << "Can't find " << filename << " at " << path << "! Program will now close." << std::endl;
//     exit(-1);
//   }
// }

//  BOOST_AUTO_TEST_CASE(  Test_ChiSqrt_Corrected )
//  {
//    // std::vector<std::vector<int>> geno;  
//    // std::vector<int> pheno;
// //   // checkInputFiles("in/fltm_CAST_0.500_imputedData.out", "geno");
// //   // checkInputFiles("in/pheno.csv", "pheno");
  
// //   // loadDataTable ( geno, "in/fltm_CAST_0.500_imputedData.out" );
// //   // loadPhenotype( pheno, "in/pheno.csv" );

//  stats::StatisticTest<stats::CHISQ> chisq;
//  for ( int i = 0; i < 20; ++i ) {
//     std::cout << i << ": " << chisq(geno[i], pheno, 3, 2) << std::endl;
//  }
//    // BOOST_CHECK_CLOSE( chisq(geno, pheno, 2, 2), 0.4233054, 0.0001 );
//  }

// BOOST_AUTO_TEST_CASE(  Test_ChiSqrt_Corrected_1 )
// {
//   std::vector<std::vector<int>> geno;  
//   // std::vector<int> pheno;
//   checkInputFiles("/home/pdt/Desktop/fltm/out/OUT/x20/fltm_CAST_0.500_imputedData.csv", "geno");
//   // checkInputFiles("in/pheno.csv", "pheno");
  
//   loadDataTable ( geno, "/home/pdt/Desktop/fltm/out/OUT/x20/fltm_CAST_0.500_imputedData.csv" );
//   // loadPhenotype( pheno, "in/pheno.csv" );
//   std::vector<int> max_vals( geno.nbrVariables(), -1);
//   // stats::StatisticTest<stats::CHISQ> chisq;
//   // for ( int i = 0; i < 20; ++i ) {
//   //    std::cout << i << ": " << chisq(geno[i], pheno, 3, 2) << std::endl;
//   // }
//   // BOOST_CHECK_CLOSE( chisq(geno, pheno, 2, 2), 0.4233054, 0.0001 );
//   std::ofstream os("in/max_vals.csv");

//   for ( int i = 0; i < geno.nbrVariables(); ++i ) {
//     for ( int j = 0; j < geno[i].nbrVariables(); ++j ) {
//       max_vals[i] = std::max( max_vals[i], geno[i][j]);
//     }
//     os << max_vals[i] << std::endl;
//   }


//   os.close();
// }




BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here
