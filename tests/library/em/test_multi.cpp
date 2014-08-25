#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE
#endif  
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <map>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <math.h>
#include <boost/lexical_cast.hpp>

#include "em/em.hpp"
#include "em/multi_em.hpp"
#include "utils/csv_parser.hpp"
#include "em/em_helper.hpp"
#include "utils/matrix_utils.hpp"
using namespace samogwas;
using namespace utility;

class Data 
{ 
};

template< typename T >
void loadGeno( std::vector< std::vector<T> >& dt,
                    const std::string& infile,
                    const char& sep = ',',
                    const char& quote = '"' ) {
  std::ifstream matrixFile(infile.c_str());
  if (!matrixFile) {
    std::cout << "not exists..." << std::endl;
    return;
  }
   dt.reserve(100000);
  utility::CSVIterator<T> matrixLine(matrixFile);
  
  for( ; matrixLine != utility::CSVIterator<T>(); ++matrixLine ) {         
    std::vector<T> row(matrixLine->size(), 0);
    for (unsigned i = 0; i < matrixLine->size(); ++i) {
      row[i] = matrixLine->at(i);
    }
    dt.push_back(row);    
  }

  dt.resize(dt.size());
}

void loadPheno( std::vector< int >& phenotype,
                    const std::string& infile ) {
  std::ifstream labPosFile(infile.c_str());
  if (!labPosFile) return;
  utility::CSVIterator<std::string> labPosLine(labPosFile);// ++labPosLine;
  for( ; labPosLine != utility::CSVIterator<std::string>(); ++labPosLine ) {    
    int pheno = boost::lexical_cast<int>( (*labPosLine)[0]);
    phenotype.push_back(pheno);
  }
}

///////////////////////////////////////////////////////////////////////////////////////
Variable createVar( const std::string lab, const int cardinality ) {
  return Variable( lab, plIntegerType(0, cardinality - 1) );
}
/////////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE( Test_NV_EM, Data ) 

BOOST_AUTO_TEST_CASE( Test_EM_Init_Values ) {
  std::vector< std::vector<int> > geno;  std::vector<int> pheno;

  loadGeno( geno, "data/em_geno_dump.csv" );
  loadPheno( pheno, "data/em_pheno_dump.csv" );
  unsigned rows = utility::nrows(geno);
  
  NV_EM em(10,1);
  ResultEM result;
  Variable Y = createVar("Y",2);
  Variable X1 = createVar("X1",2), X2 = createVar("X2", 3), X3= createVar("X3", 2);
  Variables X = X1 ^ X2 ^ X3;
  std::vector< std::vector<bool> > defTable;

  // em.run( result, Y, X, geno, defTable, 0.0000000000001);
  // unsigned N = em.theta.size(), K = em.pYX.size(), P = em.pYX[0].size();

  // BOOST_CHECK_EQUAL( N, rows );  BOOST_CHECK_EQUAL( K, 2 );  BOOST_CHECK_EQUAL( P, 3 );

  // double sum_pY = 0;
  // for ( auto& y: em.pY ) {
  //   sum_pY += y;
  //   BOOST_CHECK( y >= 0.0 );
  // }
  // // BOOST_CHECK_EQUAL( sum_pY, 1.0 );
  // BOOST_CHECK_CLOSE( sum_pY, 1.0, 0.001 );

  // for ( size_t i = 0; i < N; ++i ) {
  //   double sum_pY = 0;
  //   for ( auto& y: em.theta[i] ) {
  //     sum_pY += y;
  //     BOOST_CHECK( y >= 0.0 );
  //   }
  //   // BOOST_CHECK_EQUAL( sum_pY, 1.0 );
  //   BOOST_CHECK_CLOSE( sum_pY, 1.0, 0.001 );
  // }

  // for ( size_t y = 0; y < K; ++y ) {
  //   for ( size_t i = 0; i < P; ++i ) {
  //     double sum_pX = 0;
  //     for ( auto& x: em.pYX[y][i] )
  //     {
  //       BOOST_CHECK( x >= 0.0 );
  //       sum_pX += x;
  //     }
  //     BOOST_CHECK_CLOSE( sum_pX, 1.0, 0.001 );
  //   }
  // }
}

//////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE( Test_KL ) {
  std::vector<double> P { 0.3, 0.3, 0.4 };
  std::vector<double> Q { 0.3, 0.3, 0.4 };
  BOOST_CHECK_EQUAL( KL(P,Q), 0.0 );

  std::vector<double> Q1 { 0.4, 0.2, 0.4 };
  BOOST_CHECK( KL(P,Q1) > 0.0 );

}

void print_tabs( const std::vector< std::vector<double> >& theta,
                 const std::vector<double>& pY,
                 const std::vector< std::vector< std::vector<double> > >& pYX )
{
  for ( int i = 0; i < pY.size(); ++i ) printf("pY[%d]: %f - ", i, pY[i]);
  std::cout << std::endl;
  for ( int y = 0; y < pYX.size(); ++y ) {
    for ( int p = 0; p < pYX[y].size(); ++p ) {
      for ( int x = 0; x < pYX[y][p].size(); ++x ) {
        printf("p[%d][%d][%d]: %f, ", y, p, x, pYX[y][p][x]);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( Test_KL_Exp_Est ) {
  std::vector<double> Y1 { 0.7, 0.3 };
  std::vector< std::vector< std::vector<double> > > YX1 { { {0.3, 0.7}, {0.4, 0.3, 0.3}, {0.60, 0.40} },
                                                          { {0.7, 0.3}, {0.2, 0.2, 0.6}, {0.75, 0.25} } };
                                                        

  std::vector<double> Y2 { 0.3, 0.7 };
  std::vector< std::vector< std::vector<double> > > YX2 { { {0.7, 0.3}, {0.2, 0.2, 0.6}, {0.75, 0.25} },
                                                          { {0.3, 0.7}, {0.4, 0.3, 0.3}, {0.60, 0.40} } };
  std::vector< std::vector<int> > geno;  std::vector<int> pheno;
  loadGeno( geno, "data/em_geno_dump.csv" );
  loadPheno( pheno, "data/em_pheno_dump.csv" );
  unsigned rows = utility::nrows(geno);
  
  NV_EM em(10,1);
  ResultEM result;
  Variable Y = createVar("Y",2);
  Variable X1 = createVar("X1",2), X2 = createVar("X2", 3), X3= createVar("X3", 2);
  Variables X = X1 ^ X2 ^ X3;
  std::vector< std::vector<bool> > defTable;

  em.run( result, Y, X, geno, defTable, 0.0000000001);


  std::cout << "KL Y1:\n";
  double kl_y1 = KL( Y1, em.pY );
  double kl_yx1 = 0.0;
  for ( size_t y = 0; y < em.pYX.size(); ++y ) {
    for ( size_t p = 0; p < em.pYX[p].size(); ++p ) {
      // std::cout << "KL YX1:\n";
      printf("KL-YX1: %d,%d: %d, %d\n", y,p, (int)YX1[y][p].size(), (int)em.pYX[y][p].size());
      kl_yx1 += KL( YX1[y][p], em.pYX[y][p] );
    }
  }

  std::cout << "KL Y1:\n";
  double kl_y2 = KL( Y1, em.pY );
  double kl_yx2 = 0.0;
  for ( size_t y = 0; y < em.pYX.size(); ++y ) {
    for ( size_t p = 0; p < em.pYX[p].size(); ++p ) {
      printf("KL-YX2: %d,%d: %d, %d\n", y,p, (int)YX2[y][p].size(), (int)em.pYX[y][p].size());
      kl_yx2 += KL( YX2[y][p], em.pYX[y][p] );
    }
  }

  print_tabs( em.theta, em.pY, em.pYX );

  printf("KL1: y: %f, yx: %f\n", kl_y1, kl_yx1);
  printf("KL2: y: %f, yx: %f\n", kl_y2, kl_yx2);
  // print_tabs( em.theta, em.pY, em.pYX );



}


BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here
