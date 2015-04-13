#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE Main
#endif
#include <boost/test/unit_test.hpp>
 

#include "statistics/mutual_information.hpp" // Entropy, JointEntrop
#include <boost/algorithm/string/predicate.hpp>

#include <stdlib.h>
#include <pl.h>

using namespace samogwas;
#include <math.h> //log 
#define BOOST_TEST_MODULE Fixtures
  
class Data   
{
}; 

BOOST_FIXTURE_TEST_SUITE(Test_Entropy_Dist, Data)

 
BOOST_AUTO_TEST_CASE(EntropyCaseUniform)
{
  std::vector<int> v;
  v.push_back(0);
  v.push_back(1);
  v.push_back(2);
  v.push_back(3);

  std::vector<double> v_dist { 0.25, 0.25, 0.25, 0.25 };

  Entropy<0> entropy;
  Entropy<1> entropy_dist;
  BOOST_CHECK_EQUAL(entropy(v), log(4.0));
  BOOST_CHECK_EQUAL(entropy_dist(v_dist), log(4.0));

}


BOOST_AUTO_TEST_CASE(JointEntropy_Case_1)
{
   
  std::vector<int> varA {1,1,1,2,3,3,5}; std::vector<double> varA_dist {3.0/7, 1.0/7, 2.0/7, 1.0/7 };
  std::vector<int> varB {2,2,3,3,4,4,6}; std::vector<double> varB_dist {3.0/7, 1.0/7, 2.0/7, 1.0/7 };  
  
  Entropy<EMP> entropy; Entropy<EXACT> entropy_dist;
  double eA = entropy(varA), eA_dist = entropy_dist(varA_dist);
  double eB = entropy(varB) , eB_dist = entropy_dist(varB_dist);
  BOOST_CHECK_CLOSE( eA, 1.277034, 0.0001); BOOST_CHECK_CLOSE( eA_dist, 1.277034, 0.0001); 

  BOOST_CHECK_CLOSE( eB, 1.351784, 0.0001); BOOST_CHECK_CLOSE( eB_dist, 1.277034, 0.0001); 
  
  JointEntropy<EMP> mutEntropy; 
  double eAB = mutEntropy(varA, varB);

  BOOST_CHECK_CLOSE( eAB, 1.549826, 0.0001);

  MutualInformation<EMP> mutInfo;
  double miAB = mutInfo(varA, varB);
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

//============================================================

BOOST_AUTO_TEST_CASE(Test_Compute_Entropy_From_Dist)
{
  Entropy<0> entropy;
  std::vector<double> dist_1 { 0.5, 0.5 };
  // check
  double expected_1 = log2(2);
  BOOST_CHECK_EQUAL( entropy.compute_from_distribution(dist_1), expected_1 );

  std::vector<double> dist_3 { 0.2, 0.2, 0.2, 0.4 };  // check

  std::vector<double> dist_4 { 2, 2, 2, 4 };
  // check

  BOOST_CHECK_EQUAL(entropy.compute_from_distribution(dist_3), entropy.compute_from_distribution(dist_4));
}

/////////////////////////////////////

struct JointDistFunc {

  JointDistFunc(std::vector<std::vector<double>>& jointT): jointTab(jointT) {}
  double operator()(const plSymbol varA, const plSymbol varB,
                    const int a, const int b) const {
    return jointTab[a][b];
  }

  std::vector<std::vector<double>>& jointTab;
};

///////////////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(Test_Compute_JointEntropy_From_Dist)
{
  JointEntropy<0> jointEntropy;

  plSymbol varA("A", plIntegerType(0,2)), varB("B", plIntegerType(0,1));
  std::vector<std::vector<double>> jointTab { {1.0/6, 1.0/6}, {2.0/6, 1.0/12}, {1.0/12, 1.0/6} };
  JointDistFunc jf(jointTab);
  double val = jointEntropy.compute_from_joint_function(varA,varB,jf);

  Entropy<0> entropy;
  std::vector<double> dist { 1.0/6, 1.0/6, 2.0/6, 1.0/12, 1.0/12, 1.0/6 };
  double expected = entropy.compute_from_distribution(dist);

  BOOST_CHECK_EQUAL( val, expected );
}



// w

// bool JointFunc::is_latent(const plSymbol varA) {
//   return boost::starts_with(varA.name(), "latent");
// }



BOOST_AUTO_TEST_SUITE_END()  /// Test Info_Theo ends here
// ///////////////////////////////////////////////////////////////////////////////////










