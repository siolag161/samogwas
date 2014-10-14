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
#include "gwas/gwas_basic_strategy.hpp"

using namespace stats;
class Data
{
};
BOOST_FIXTURE_TEST_SUITE( Test_Basic_Strategy, Data )


//////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE( Test_Basic ) {

    }


BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here





