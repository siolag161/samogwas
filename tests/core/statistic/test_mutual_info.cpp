#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE Main
#endif
#include <boost/test/unit_test.hpp>
 

#include "statistics/mutual_information.hpp" // Entropy, JointEntropy

#include <stdlib.h>

using namespace samogwas;
#include <math.h> //log 
#define BOOST_TEST_MODULE Fixtures
  
class Data   
{
}; 

BOOST_FIXTURE_TEST_SUITE(Test_Entropy, Data)

 
BOOST_AUTO_TEST_CASE(EntropyCaseUniform)
{
  std::vector<int> v;
  v.push_back(1);
  v.push_back(2);
  v.push_back(3);
  v.push_back(4);
  v.push_back(-1);

  Entropy<0>entropy;
  BOOST_CHECK_EQUAL(entropy(v), log(4.0));
}


BOOST_AUTO_TEST_CASE(JointEntropy_Case_1)
{
   
  std::vector<int> varA;
  varA.push_back(1);
  varA.push_back(1);
  varA.push_back(1);
  varA.push_back(2);

  varA.push_back(-1);

  varA.push_back(3);
  varA.push_back(3);
  varA.push_back(5);

  std::vector<int> varB;
  varB.push_back(2);
  varB.push_back(2);
  varB.push_back(3);
  varB.push_back(3);

  varB.push_back(-1);
  
  varB.push_back(4);
  varB.push_back(4);
  varB.push_back(6);
  
  
  Entropy<EMP> entropy;
  double eA = entropy(varA);
  double eB = entropy(varB);
  BOOST_CHECK_CLOSE( eA, 1.277034, 0.0001);
  BOOST_CHECK_CLOSE( eB, 1.351784, 0.0001);
  
  JointEntropy<EMP> mutEntropy;
  double eAB = mutEntropy(varA, varB);
  BOOST_CHECK_CLOSE( eAB, 1.549826, 0.0001);

  MutualInformation<EMP> mutInfo;
  double miAB = mutInfo(varA, varB);
  // BOOST_CHECK_EQUAL((eA + eB - eAB), miAB);
  BOOST_CHECK(std::abs(eA + eB - eAB - miAB) < 0.0001); 
}

BOOST_AUTO_TEST_CASE(JointEntropy_Case_2)
{
   
  std::vector<int> varA;
  varA.push_back(1);
  varA.push_back(2);
  varA.push_back(2);
  varA.push_back(3);
  varA.push_back(1);
  varA.push_back(3);

  std::vector<int> varB;
  varB.push_back(2);
  varB.push_back(2);
  varB.push_back(3);
  varB.push_back(3);
  varB.push_back(2);
  varB.push_back(3);

  Entropy<EMP> entropy;
  double eA = entropy(varA);
  double eB = entropy(varB);
  BOOST_CHECK_CLOSE( eA, 1.098612, 0.0001);
  BOOST_CHECK_CLOSE( eB, 0.6931472, 0.0001);
  
  JointEntropy<EMP> mutEntropy;
  double eAB = mutEntropy(varA, varB);
  BOOST_CHECK_CLOSE( eAB, 1.329661, 0.0001);

  MutualInformation<EMP> mutInfo;
  double miAB = mutInfo(varA, varB);
  // BOOST_CHECK_EQUAL((eA + eB - eAB), miAB);
 
}


BOOST_AUTO_TEST_CASE(VectorVectorMutInfo_Iterator)
{
  std::vector<int> varA;
  varA.push_back(1);
  varA.push_back(2);
  varA.push_back(3);
  varA.push_back(1);
  varA.push_back(5);
  varA.push_back(6);

  std::vector<int> varB;
  varB.push_back(1);
  varB.push_back(11);
  varB.push_back(12);
  varB.push_back(10);
  varB.push_back(11);
  varB.push_back(13);
   
  MutualInformation<EMP> mutInfo;  
  BOOST_CHECK_CLOSE( mutInfo(varA.begin(), varA.end(), varB.begin()), 1.329661, 0.0001);
 
}


BOOST_AUTO_TEST_CASE(VectorVectorMutInfo_stdVec)
{ 
  std::vector<int> varB;
  varB.push_back(1); varB.push_back(2); varB.push_back(3);
  varB.push_back(1); varB.push_back(5); varB.push_back(6);
  
  std::vector<int> varA;
  varA.push_back(1);varA.push_back(11); varA.push_back(12);
  varA.push_back(10);varA.push_back(11); varA.push_back(13);    
 
  MutualInformation<EMP> mutInfo;
  
  BOOST_CHECK_CLOSE( mutInfo(varA, varB), 1.329661,  0.0001);
  BOOST_CHECK_EQUAL(mutInfo(varA, varB), mutInfo(varA, varB));
 
}


///////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(MutInfo_stdVec)
{
  
 std::vector<int> varA;
  varA.push_back(1);
  varA.push_back(2);
  varA.push_back(3);
  varA.push_back(1);
  varA.push_back(5);
  varA.push_back(6);

  std::vector<int> varB;
  varB.push_back(1);
  varB.push_back(11);
  varB.push_back(12);
  varB.push_back(10);
  varB.push_back(11);
  varB.push_back(13);  
  
  std::vector< std::vector<int> >  mat;
  mat.push_back(varA);
  mat.push_back(varB);    
     
  MutualInformation<EMP> mutInfo; 

}



BOOST_AUTO_TEST_SUITE_END()  /// Test Info_Theo ends here
// ///////////////////////////////////////////////////////////////////////////////////
