#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE
#endif  
#include <boost/test/unit_test.hpp>

#include <omp.h>

#include <vector>
#include "Distance.hpp"
#include "chrono_measure.hpp"

using namespace clustering; 
class Data 
{ 

};
 
BOOST_FIXTURE_TEST_SUITE( Test_Para_Entropy, Data ) 
 
BOOST_AUTO_TEST_CASE( Test_Entropy ) {


}

BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here
