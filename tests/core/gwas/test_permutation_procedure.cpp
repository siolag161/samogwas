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

#include "statistics/permutation_test.hpp"

using namespace stats;
class Data
{
};
BOOST_FIXTURE_TEST_SUITE( Test_GWAS, Data )


//////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE( Test_Def ) {
  std::vector<double> distri;
  std::vector<double> pvals;
  StatTest* statTest;


}


BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here





