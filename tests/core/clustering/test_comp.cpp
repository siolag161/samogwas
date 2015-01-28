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
#include <memory>

#include "clustering/compare_measure.hpp"
#include "clustering/rand.hpp"
#include "clustering/mirkin.hpp"
#include "clustering/clustering_quality.hpp"

#include "test_clustering.hpp"

class Data 
{ 
};

BOOST_FIXTURE_TEST_SUITE( Test_Clustering_Comp, Data ) 

BOOST_AUTO_TEST_CASE( Test_Rand ) {
  Clustering cl1, cl2;
  std::vector<int> c11 { 0,1,2 };
  std::vector<int> c12 { 3,4,5 };
  std::vector<int> c13 { 6,7,8 };

  std::vector<int> c21 { 1,3,4,5,7,8 };
  std::vector<int> c22 { 0,2,6 };

  cl1.push_back(c11); cl1.push_back(c12); cl1.push_back(c13);
  cl2.push_back(c21); cl2.push_back(c22); // cl2.push_back(c23);

  samogwas::RandIndex rand;
  double rv = rand(cl1, cl2);
  BOOST_CHECK_CLOSE( rv, 0.5277778, 0.0001 );
}


BOOST_AUTO_TEST_CASE( Test_Rand_Partition ) {
  // Clustering cl1, cl2;
  // true.id <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  // pred.id <- c(2, 1, 2, 1, 1, 1, 2, 1, 1)
  // label   <- c(0, 0, 0, 0, 1, 0, 2, 0, 0)
  Partition pA, pB;
  pA.setLabel(0,0).setLabel(1,0).setLabel(2,0).setLabel(3,1).setLabel(4,1).setLabel(5,1).setLabel(6,2).setLabel(7,2).setLabel(8,2);
  pB.setLabel(0,1).setLabel(1,0).setLabel(2,1).setLabel(3,0).setLabel(4,0).setLabel(5,0).setLabel(6,1).setLabel(7,0).setLabel(8,0);  

  // cl1.push_back(c11); cl1.push_back(c12); cl1.push_back(c13);
  // cl2.push_back(c21); cl2.push_back(c22); // cl2.push_back(c23);

  samogwas::RandIndex r;
  double rv = r(pA,pB);
  BOOST_CHECK_CLOSE( rv, 0.5277778, 0.0001 );

}

BOOST_AUTO_TEST_CASE( Test_Adj_Rand_Partition ) {
  Partition pA, pB;
  pA.setLabel(0,0).setLabel(1,0).setLabel(2,0).setLabel(3,1).setLabel(4,1).setLabel(5,1).setLabel(6,2).setLabel(7,2).setLabel(8,2);
  pB.setLabel(0,1).setLabel(1,0).setLabel(2,1).setLabel(3,0).setLabel(4,0).setLabel(5,0).setLabel(6,1).setLabel(7,0).setLabel(8,0);  

  samogwas::AdjustedRandIndex arand;
  double rv = arand(pA,pB);
  BOOST_CHECK_CLOSE( rv, 0.05555556, 0.0001 ); 
}

BOOST_AUTO_TEST_CASE( Test_Mirkin ) {
  // inline double mirkin_distance( const Clustering& c1, const Clustering& c2 ) {
  samogwas::Clustering cl1, cl2, cl3;

  std::vector<int> c1 {1,3,5,7};
  std::vector<int> c2 {0,2,4,6};

  std::vector<int> c3 {1,2,5,7};
  std::vector<int> c4 {0,3,4,6};

  cl1.push_back(c1); cl1.push_back(c2);
  cl2.push_back(c1); cl2.push_back(c2);
  cl3.push_back(c3); cl3.push_back(c4);
  
  samogwas::Mirkin mirkin;
  double mv = mirkin(cl1, cl2);
  BOOST_CHECK_EQUAL( mv, 0.0 );
  BOOST_CHECK( mirkin(cl2,cl3) != 0 );

}

BOOST_AUTO_TEST_CASE( Test_Clustering_Quality ) {
  samogwas::Clustering cl1, cl2;
  std::vector<int> c11 { 0,1,2 };
  std::vector<int> c12 { 3,4,5 };
  std::vector<int> c13 { 6,7,8 };

  std::vector<int> c21 { 1,3,4,5,7,8 };
  std::vector<int> c22 { 0,2,6 };

  cl1.push_back(c11); cl1.push_back(c12); cl1.push_back(c13);
  cl2.push_back(c21); cl2.push_back(c22); // cl2.push_back(c23);

  // samogwas::RandIndex rand;
  // double rv = rand(cl1, cl2);
  // BOOST_CHECK_CLOSE( rv, 0.5277778, 0.0001 );

  samogwas::SpatialCoherenceCriteria spatial;

  samogwas::Partition p1(cl1);
  samogwas::Partition p2(cl2);

  BOOST_CHECK_EQUAL( spatial.nbrLabelChanges(p1), 2.0);
  BOOST_CHECK_EQUAL( spatial.nbrLabelChanges(p2), 5.0);

  BOOST_CHECK_EQUAL( spatial.expectedNbrLabelChanges(cl1), 6.0);
  BOOST_CHECK_EQUAL( spatial.expectedNbrLabelChanges(cl2), 4.0);

  BOOST_CHECK_EQUAL( spatial.compute(cl1), spatial.compute(p1) );
  BOOST_CHECK( spatial.compute(cl1) > spatial.compute(cl2) );


}


BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here
