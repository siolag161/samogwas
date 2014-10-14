/*
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
#include "em/naive_bayes_em.hpp"
#include "utils/csv_parser.hpp"
#include "em/em_helper.hpp"
#include "utils/matrix_utils.hpp"
#include "fltm/fltm.hpp"
#include "data_generation.hpp"

using namespace samogwas;
using namespace utility;

class Data 
{ 
};

////////////////////////////////////////////////
BOOST_FIXTURE_TEST_SUITE( Test_Naive_Bayes_EM, Data ) 
/////////////////////////////////////////////////

template< typename T >
void loadGeno( std::vector< std::vector<T> >& dt,
                    const std::string& infile,
                    const char& sep = ',',
                    const char& quote = '"' ) {
  std::ifstream matrixFile(infile.c_str());
  if (!matrixFile) {
    std::cout << infile <<": not exists..." << std::endl;
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


BOOST_AUTO_TEST_CASE( Test_EM_Init_Values ) {
  std::vector< std::vector<int> > geno;  std::vector<int> pheno;

  loadGeno( geno, "data/em_geno_dump.csv" );
  loadPheno( pheno, "data/em_pheno_dump.csv" );
  unsigned rows = utility::nrows(geno);

  size_t nclusts = 5, ncols = 40;
  size_t N = 3, CARD = 3, MAX_POS = 50;
  int nrows = nclusts*N;
  std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  auto data = data_gen::GenerateClusteredData( nclusts, N, CARD, ncols )();

        EMInterface* em = new NaiveBayesEM(10,1);
  ResultEM result;
  Variable Y = createVar("Y", CARD);
  Variable X1 = createVar("X1", CARD), X2 = createVar("X2", CARD), X3= createVar("X3", CARD);
  Variables X = X1 ^ X2 ^ X3;
  std::vector< std::vector<bool> > defTable;


  em->run(result, Y, X, data, 0.01);

}


BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here
*/
